"""Run a reproducible grid of individual POSYDON formation experiments."""

import argparse
import copy
import hashlib
import importlib.metadata
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from formation_diagnostics import (
    BINARY_FIELDS, STAR_FIELDS, json_safe, make_diagnostic_ce_class, number,
    orbital_screen, snapshot, summarize_history,
)


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_STUDY = PROJECT_DIR / "config" / "targeted_formation.json"
CE_OVERRIDE_KEYS = {
    "common_envelope_efficiency", "common_envelope_option_for_lambda",
    "common_envelope_lambda_default", "common_envelope_option_for_HG_star",
    "common_envelope_alpha_thermal", "core_definition_H_fraction",
    "core_definition_He_fraction", "common_envelope_option_after_succ_CEE",
}


def load_study(args):
    spec = json.loads(args.study.read_text())
    for key in ("primary_masses", "companion_masses", "periods_days", "seed",
                "kick_directions"):
        value = getattr(args, key, None)
        if value is not None:
            spec[key] = value
    if args.no_kick_replays:
        spec["kick_speeds_kms"] = []
    if args.scenario:
        known = {s["name"] for s in spec["scenarios"]}
        unknown = set(args.scenario) - known
        if unknown:
            raise ValueError(f"Unknown scenarios: {sorted(unknown)}")
        spec["scenarios"] = [s for s in spec["scenarios"] if s["name"] in args.scenario]
    for key in ("primary_masses", "companion_masses", "periods_days"):
        values = spec[key]
        if not values or len(set(values)) != len(values):
            raise ValueError(f"{key} must contain distinct values")
        if not all(math.isfinite(v) and v > 0 for v in values):
            raise ValueError(f"{key} must contain positive finite values")
    if max(spec["companion_masses"]) > min(spec["primary_masses"]):
        raise ValueError("Keep companion masses <= all primary masses in this grid")
    if not 0 <= spec["eccentricity"] < 1:
        raise ValueError("eccentricity must be in [0, 1)")
    if not math.isfinite(spec["max_time_yr"]) or spec["max_time_yr"] <= 0:
        raise ValueError("max_time_yr must be positive and finite")
    if not math.isfinite(spec["metallicity_zsun"]) or spec["metallicity_zsun"] <= 0:
        raise ValueError("metallicity_zsun must be positive and finite")
    if not isinstance(spec["seed"], int) or spec["seed"] < 0:
        raise ValueError("seed must be a nonnegative integer")
    if not isinstance(spec["kick_directions"], int) or spec["kick_directions"] < 1:
        raise ValueError("kick_directions must be a positive integer")
    if not all(math.isfinite(v) and v >= 0 for v in spec["kick_speeds_kms"]):
        raise ValueError("kick_speeds_kms must contain nonnegative finite values")
    names = []
    for scenario in spec["scenarios"]:
        name = scenario["name"]
        if not re.fullmatch(r"[a-zA-Z0-9_-]+", name):
            raise ValueError(f"Invalid scenario name: {name}")
        names.append(name)
        overrides = scenario["ce_overrides"]
        if set(overrides) - CE_OVERRIDE_KEYS:
            raise ValueError(f"Unsupported CE override keys in {name}")
        alpha = overrides.get("common_envelope_efficiency", 1.0)
        if not math.isfinite(alpha) or alpha <= 0:
            raise ValueError("CE efficiency must be positive and finite")
        if overrides.get("core_definition_H_fraction", 0.3) not in (0.01, 0.1, 0.3):
            raise ValueError("POSYDON supports H core boundaries 0.01, 0.1, or 0.3")
    if not names or len(names) != len(set(names)):
        raise ValueError("Scenarios must have distinct names")
    screen = spec["orbital_screen"]
    for low, high in (("companion_mass_min", "companion_mass_max"),
                      ("period_min_days", "period_max_days")):
        if not (math.isfinite(screen[low]) and math.isfinite(screen[high])
                and 0 < screen[low] <= screen[high]):
            raise ValueError(f"Invalid orbital screen: {low}, {high}")
    population_config = args.population_config or PROJECT_DIR / spec["population_config"]
    population_config = population_config.resolve()
    if not population_config.is_file():
        raise FileNotFoundError(population_config)
    spec["population_config"] = str(population_config)
    return spec


def initial_cases(spec):
    return [{"case_id": i, "M1_initial": m1, "M2_initial": m2,
             "period_initial_days": period, "q_initial": m2 / m1}
            for i, (m1, m2, period) in enumerate(itertools.product(
                spec["primary_masses"], spec["companion_masses"], spec["periods_days"]))]


def support_reason(case, limits):
    if not limits["m1_min"] <= case["M1_initial"] <= limits["m1_max"]:
        return "primary_mass_outside_HMS_HMS_grid"
    minimum_q = max(limits["q_min"], 0.5 / case["M1_initial"])
    if not minimum_q <= case["q_initial"] <= limits["q_max"]:
        return "mass_ratio_outside_HMS_HMS_grid"
    if not limits["p_min"] <= case["period_initial_days"] <= limits["p_max"]:
        return "period_outside_HMS_HMS_grid"
    return None


def grid_limits_from_file(spec):
    """Fast preflight of the ZAMS grid without loading all interpolators."""
    import ast
    from configparser import ConfigParser
    import h5py
    from posydon.config import PATH_TO_POSYDON_DATA
    from posydon.utils.common_functions import convert_metallicity_to_string

    config = ConfigParser()
    config.read(spec["population_config"])
    grid_directory = ast.literal_eval(config.get("grid_paths", "HMS_HMS_path", fallback="None"))
    grid_directory = Path(grid_directory or Path(PATH_TO_POSYDON_DATA) / "HMS-HMS")
    grid_file = grid_directory / f"{convert_metallicity_to_string(spec['metallicity_zsun'])}_Zsun.h5"
    with h5py.File(grid_file, "r") as grid:
        values = grid["grid/initial_values"][:]
    # These mass-ratio bounds match MS_MS_step in the installed POSYDON 2.3.4.
    return {"m1_min": float(values["star_1_mass"].min()),
            "m1_max": float(values["star_1_mass"].max()), "q_min": 0.05, "q_max": 1.0,
            "p_min": float(values["period_days"].min()),
            "p_max": float(values["period_days"].max())}


def write_json(path, value):
    path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + "\n")


def append_csv(path, records):
    if records:
        frame = pd.DataFrame(records)
        frame.reindex(columns=sorted(frame.columns)).to_csv(
            path, mode="a", header=not path.exists(), index=False)


def make_properties(spec):
    from posydon.binary_evol.simulationproperties import SimulationProperties
    from posydon.popsyn.io import simprop_kwargs_from_ini

    kwargs = simprop_kwargs_from_ini(spec["population_config"])
    kwargs["step_CE"] = (make_diagnostic_ce_class(), kwargs["step_CE"][1])
    properties = SimulationProperties(**kwargs)
    properties.max_simulation_time = spec["max_time_yr"]
    properties.end_events = ["CC1", "CC2"]
    properties.end_states = ["merged", "disrupted", "initial_RLOF"]
    properties.load_steps(metallicity=spec["metallicity_zsun"],
                          RNG=np.random.default_rng(spec["seed"]))
    if not properties.step_SN.kick:
        properties.close()
        raise ValueError("Use step_SN kick=True for impulsive mass-loss dynamics; "
                         "the zero-kick control prescribes a velocity of zero.")
    return properties


def make_binary(case, properties, spec):
    from posydon.binary_evol.binarystar import BinaryStar
    from posydon.binary_evol.singlestar import SingleStar

    return BinaryStar(
        SingleStar(mass=case["M1_initial"], state="H-rich_Core_H_burning"),
        SingleStar(mass=case["M2_initial"], state="H-rich_Core_H_burning"),
        index=case["case_id"], properties=properties,
        time=0.0, state="detached", event="ZAMS",
        orbital_period=case["period_initial_days"], eccentricity=spec["eccentricity"])


def history_frame(binary):
    columns = list(BINARY_FIELDS) + ["step_names"]
    columns += [f"S{i}_{key}" for i in (1, 2) for key in STAR_FIELDS]
    return binary.to_df(extra_columns={"step_names": "string"})[columns].reset_index(drop=True)


def kick_replays(pre_sn, case, scenario, spec):
    """Reuse a fixed pre-collapse binary; prescribed speeds bypass kick scaling."""
    records = []
    component = 1 if pre_sn.event == "CC1" else 2
    rng = np.random.default_rng(np.random.SeedSequence([spec["seed"], case["case_id"], 1]))
    for speed in spec["kick_speeds_kms"]:
        trials = (1 if speed == 0 and pre_sn.eccentricity == 0
                  else spec["kick_directions"])
        for direction in range(trials):
            trial = copy.deepcopy(pre_sn, {id(pre_sn.properties): pre_sn.properties})
            star = getattr(trial, f"star_{component}")
            phi = rng.uniform(0, 2 * np.pi)
            theta = np.arccos(rng.uniform(-1, 1))
            mean_anomaly = rng.uniform(0, 2 * np.pi)
            star.natal_kick_velocity = float(speed)
            star.natal_kick_azimuthal_angle = float(phi)
            star.natal_kick_polar_angle = float(theta)
            star.natal_kick_mean_anomaly = float(mean_anomaly)
            record = {"scenario": scenario, "case_id": case["case_id"],
                      "kick_speed_kms": speed, "direction_id": direction,
                      "azimuth_rad": phi, "polar_rad": theta,
                      "mean_anomaly_rad": mean_anomaly, "error": ""}
            try:
                # Keep kick=True: disabling the kick step also changes its mass-loss dynamics.
                trial.properties.step_SN(trial)
                record.update(orbital_screen(trial, spec["orbital_screen"]))
                record.update(state=trial.state, period_days=number(trial.orbital_period),
                              eccentricity=number(trial.eccentricity))
            except Exception as exc:
                record.update(error=f"{type(exc).__name__}: {exc}", state="error",
                              bound_ns_ms=False, orbital_screen_pass=False,
                              periastron_detached=False, M_NS=None, M_companion=None,
                              period_days=None, eccentricity=None)
            records.append(record)
    return records


def run_case(case, scenario, properties, spec, output, limits):
    started = time.monotonic()
    result = {**case, "scenario": scenario, "status": "", "error": "",
              "ce_count": 0, "ce_survivals": 0, "reached_first_collapse": False,
              "ce_early_low_mass_attempts": 0, "ce_early_low_mass_survivals": 0,
              "bound_ns_ms": False, "orbital_screen_pass": False,
              "early_stripped_orbital_screen_pass": False}
    reason = support_reason(case, limits)
    if reason:
        result.update(status="unsupported_initial_conditions", error=reason,
                      bottleneck="initial_grid_coverage", elapsed_seconds=0.0)
        return result
    binary = make_binary(case, properties, spec)
    binary.study_ce_records = []
    properties.step_SN.RNG = np.random.default_rng(np.random.SeedSequence(
        [spec["seed"], case["case_id"], 0]))
    pre_sn = None
    try:
        binary.evolve()
        if binary.event in ("CC1", "CC2"):
            result["reached_first_collapse"] = True
            component = 1 if binary.event == "CC1" else 2
            result["collapse_component"] = component
            result.update({f"preSN_{key}": value for key, value in snapshot(binary).items()})
            pre_sn = copy.deepcopy(binary, {id(properties): properties})
            binary.run_step()
            result["remnant_state"] = getattr(binary, f"star_{component}").state
            result["SN_type"] = getattr(getattr(binary, f"star_{component}"), "SN_type", None)
            result.update(orbital_screen(binary, spec["orbital_screen"]))
            result["status"] = "first_collapse_processed"
        else:
            result["status"] = "terminated_before_first_collapse"
    except Exception as exc:
        result.update(status="simulation_error", error=f"{type(exc).__name__}: {exc}")
        (output / "errors" / f"{scenario}_{case['case_id']:04d}.txt").write_text(traceback.format_exc())
    finally:
        history = history_frame(binary)
        history.to_csv(output / "histories" / f"{scenario}_{case['case_id']:04d}.csv", index=False)
        result.update(summarize_history(history))
        result.update(stop_time_yr=binary.time, stop_state=binary.state,
                      stop_event=binary.event, stop_period_days=number(binary.orbital_period),
                      stop_eccentricity=number(binary.eccentricity))
        ce_records = binary.study_ce_records
        result["ce_count"] = len(ce_records)
        result["ce_survivals"] = sum(r["after"]["state"] == "detached"
                                     and "error" not in r for r in ce_records)
        for record in ce_records:
            donor = record["donor"]
            companion = 3 - donor
            before = record["before"]
            early_low_mass = (
                number(before[f"S{donor}_center_h1"]) <= 0.01
                and number(before[f"S{donor}_center_he4"]) >= 0.9
                and before[f"S{companion}_state"] == "H-rich_Core_H_burning"
                and spec["orbital_screen"]["companion_mass_min"]
                    <= number(before[f"S{companion}_mass"])
                    <= spec["orbital_screen"]["companion_mass_max"])
            result["ce_early_low_mass_attempts"] += int(early_low_mass)
            result["ce_early_low_mass_survivals"] += int(
                early_low_mass and record["after"]["state"] == "detached"
                and "error" not in record)
        if ce_records:
            first = ce_records[0]
            donor = first["donor"]
            result.update(first_ce_time_yr=first["before"]["time"],
                          first_ce_donor=donor,
                          first_ce_donor_state=first["before"][f"S{donor}_state"],
                          first_ce_center_he4=first["before"][f"S{donor}_center_he4"])
        with (output / "ce_events.jsonl").open("a") as stream:
            for event_index, record in enumerate(ce_records):
                stream.write(json.dumps(json_safe({"scenario": scenario,
                    "case_id": case["case_id"], "ce_event": event_index, **record}),
                    allow_nan=False) + "\n")
        component = result.get("collapse_component", 1)
        early = result.get(f"S{component}_stripping_timing", "").startswith(
            ("before_helium_ignition", "early_helium"))
        result["early_stripped_orbital_screen_pass"] = bool(result["orbital_screen_pass"] and early)
    if pre_sn is not None and result.get("remnant_state") == "NS" and not result["error"]:
        append_csv(output / "kick_trials.csv", kick_replays(pre_sn, case, scenario, spec))
    result["bottleneck"] = case_bottleneck(result, spec["orbital_screen"])
    result["elapsed_seconds"] = round(time.monotonic() - started, 4)
    return result


def case_bottleneck(result, screen):
    if result["error"]:
        return "simulation_error"
    if result["stop_state"] == "merged":
        return "ce_merger" if result["ce_count"] else "other_merger"
    if not result["reached_first_collapse"]:
        return "terminated_before_collapse"
    if result.get("remnant_state") != "NS":
        return "non_NS_remnant"
    if result["stop_state"] == "disrupted":
        return "SN_disruption"
    if result["orbital_screen_pass"]:
        return "birth_orbit_screen_pass"
    if result["bound_ns_ms"]:
        if not result["periastron_detached"]:
            return "post_SN_periastron_overflow"
        if not screen["companion_mass_min"] <= result["M_companion"] <= screen["companion_mass_max"]:
            return "companion_mass_outside_screen"
        if result["stop_period_days"] < screen["period_min_days"]:
            return "post_SN_orbit_too_short"
        return "post_SN_orbit_too_wide"
    return "no_bound_NS_MS"


def save_summary(results, output, spec):
    frame = pd.DataFrame(results)
    frame.to_csv(output / "cases.csv", index=False)
    summary = {"cases": len(frame), "scenarios": {},
               "interpretation": "Unweighted controlled experiments, stopped after the first collapse. "
               "The orbital screen is a birth-orbit proxy, not present-day Gaia detectability. "
               "Early stripping support does not validate an Ertl-grid correspondence."}
    for scenario, group in frame.groupby("scenario", sort=False):
        summary["scenarios"][scenario] = {
            "cases": len(group), "statuses": group.status.value_counts().to_dict(),
            "bottlenecks": group.bottleneck.value_counts().to_dict(),
            "stop_states": group.stop_state.value_counts().to_dict() if "stop_state" in group else {},
            "ce_events": int(group.ce_count.sum()), "ce_survivals": int(group.ce_survivals.sum()),
            "ce_early_low_mass_attempts": int(group.ce_early_low_mass_attempts.sum()),
            "ce_early_low_mass_survivals": int(group.ce_early_low_mass_survivals.sum()),
            "first_collapse_records": int(group.reached_first_collapse.sum()),
            "bound_ns_ms_at_birth": int(group.bound_ns_ms.sum()),
            "orbital_screen_pass": int(group.orbital_screen_pass.sum()),
            "early_stripped_orbital_screen_pass": int(group.early_stripped_orbital_screen_pass.sum()),
        }
    kick_path = output / "kick_trials.csv"
    if kick_path.exists():
        trials = pd.read_csv(kick_path)
        summaries = []
        for (scenario, case_id, speed), group in trials.groupby(
                ["scenario", "case_id", "kick_speed_kms"], sort=False):
            valid = group[group.error.isna() | group.error.eq("")]
            summaries.append({"scenario": scenario, "case_id": case_id,
                "kick_speed_kms": speed, "trials": len(group), "errors": len(group) - len(valid),
                "bound_ns_ms_fraction": valid.bound_ns_ms.mean(),
                "orbital_screen_fraction": valid.orbital_screen_pass.mean()})
        pd.DataFrame(summaries).to_csv(output / "kick_summary.csv", index=False)
    ce_path = output / "ce_events.jsonl"
    if ce_path.exists():
        events = [json.loads(line) for line in ce_path.read_text().splitlines()]
        if events:
            pd.json_normalize(events, sep="_").to_csv(output / "ce_events.csv", index=False)
    write_json(output / "summary.json", summary)
    report = ["# Targeted formation results", "", summary["interpretation"], "",
              "| CE scenario | Cases | Unsupported | Errors | CE survivals | First collapses | Bound NS + MS at birth | Orbit screen | Early + orbit screen |",
              "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for name, counts in summary["scenarios"].items():
        statuses = counts["statuses"]
        report.append(f"| {name} | {counts['cases']} | "
                      f"{statuses.get('unsupported_initial_conditions', 0)} | "
                      f"{statuses.get('simulation_error', 0)} | {counts['ce_survivals']} | "
                      f"{counts['first_collapse_records']} | {counts['bound_ns_ms_at_birth']} | "
                      f"{counts['orbital_screen_pass']} | {counts['early_stripped_orbital_screen_pass']} |")
    report += ["", "Counts above use one baseline kick per binary. Inspect `kick_summary.csv` for",
               "replays at prescribed speeds, and `ce_events.csv` for envelope energetics.", "",
               "The orbital screen applies only the configured companion mass and period",
               "bounds and detachment at periastron. It does not apply the observed NS-mass",
               "or eccentricity cuts, present-day evolution, or Gaia detection probabilities.", "",
               "The alpha estimate is for fast-phase Roche detachment at fixed masses/radii;",
               "it is not the full CE-survival threshold."]
    (output / "report.md").write_text("\n".join(report) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=DEFAULT_STUDY)
    parser.add_argument("--population-config", type=Path)
    parser.add_argument("--primary-masses", type=float, nargs="+")
    parser.add_argument("--companion-masses", type=float, nargs="+")
    parser.add_argument("--periods-days", type=float, nargs="+")
    parser.add_argument("--scenario", action="append", help="Select a named scenario; may be repeated")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--kick-directions", type=int)
    parser.add_argument("--no-kick-replays", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ns-progenitors-matplotlib"))
    try:
        spec = load_study(args)
        limits = grid_limits_from_file(spec)
    except (ValueError, KeyError, OSError) as exc:
        parser.error(str(exc))
    cases = initial_cases(spec)
    unsupported = sum(support_reason(case, limits) is not None for case in cases)
    total = len(cases) * len(spec["scenarios"])
    print(f"{len(cases)} initial binaries x {len(spec['scenarios'])} CE scenarios = {total} cases", flush=True)
    print(f"{unsupported} initial binaries outside ZAMS grid bounds; these will be flagged, not evolved.", flush=True)
    if args.dry_run:
        print(json.dumps({"study": spec, "ZAMS_grid_limits": limits}, indent=2))
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output = (args.output or PROJECT_DIR / "outputs" / "targeted_formation" / stamp).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"Output directory is not empty: {output}")
    for directory in (output, output / "histories", output / "errors"):
        directory.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output / "matplotlib_cache"))
    shutil.copyfile(spec["population_config"], output / "population_config.ini")
    write_json(output / "manifest.json", {
        "study": spec, "ZAMS_grid_limits": limits,
        "posydon_version": importlib.metadata.version("posydon"),
        "population_config_sha256": hashlib.sha256(Path(spec["population_config"]).read_bytes()).hexdigest(),
        "code_sha256": {name: hashlib.sha256((PROJECT_DIR / name).read_bytes()).hexdigest()
                        for name in ("run_targeted.py", "formation_diagnostics.py")},
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "sampling": "unweighted Cartesian grid; no star-formation-time sampling",
        "stopping_rule": "first collapse and its immediate outcome, or earlier merger/termination",
        "alpha_estimate": "Fixed-mass/radius fast-detachment estimate; excludes subsequent "
                          "mass transfer and structural response. Not a full survival threshold.",
    })
    print(f"Outputs: {output}\nLoading POSYDON steps...", flush=True)
    properties = make_properties(spec)
    ce = properties.step_CE
    baseline = {key: getattr(ce, key) for key in CE_OVERRIDE_KEYS}
    manifest = json.loads((output / "manifest.json").read_text())
    manifest["effective_ce_settings"] = {s["name"]: {**baseline, **s["ce_overrides"]}
                                        for s in spec["scenarios"]}
    manifest["supernova_settings"] = {
        key: getattr(properties.step_SN, key) for key in
        ("mechanism", "engine", "kick", "kick_normalisation", "kick_prescription",
         "sigma_kick_CCSN_NS", "sigma_kick_ECSN")}
    write_json(output / "manifest.json", manifest)
    actual_limits = {key: float(getattr(properties.step_HMS_HMS, key)) for key in limits}
    if actual_limits != limits:
        properties.close()
        raise RuntimeError("Preflight bounds disagree with the installed ZAMS step")
    results = []
    try:
        for scenario in spec["scenarios"]:
            for key, value in {**baseline, **scenario["ce_overrides"]}.items():
                setattr(ce, key, value)
            for case in cases:
                result = run_case(case, scenario["name"], properties, spec, output, limits)
                results.append(result)
                save_summary(results, output, spec)
                print(f"[{len(results)}/{total}] {scenario['name']}: "
                      f"M1={case['M1_initial']:g}, M2={case['M2_initial']:g}, "
                      f"P={case['period_initial_days']:g} d -> "
                      f"{result['status']} / {result.get('stop_state', 'unsupported')}", flush=True)
    finally:
        properties.close()
    print(f"Completed. Read {output / 'cases.csv'} and {output / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
