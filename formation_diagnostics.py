"""Diagnostics for a controlled binary-formation study, without population weights."""

import math

import numpy as np


MS_STATE = "H-rich_Core_H_burning"
STAR_FIELDS = (
    "state", "mass", "log_R", "he_core_mass", "he_core_radius",
    "co_core_mass", "co_core_radius", "center_h1", "center_he4",
    "surface_h1", "surface_he4", "log_LH", "log_LHe",
)
BINARY_FIELDS = (
    "time", "state", "event", "separation", "orbital_period", "eccentricity",
    "mass_transfer_case",
)


def number(value):
    """Represent missing/nonfinite stellar attributes as NaN for CSV output."""
    try:
        result = float(value)
        return result if math.isfinite(result) else math.nan
    except (TypeError, ValueError):
        return math.nan


def snapshot(binary):
    result = {key: getattr(binary, key, None) for key in BINARY_FIELDS}
    for component in (1, 2):
        star = getattr(binary, f"star_{component}")
        result.update({f"S{component}_{key}": getattr(star, key, None)
                       for key in STAR_FIELDS})
    return result


def stripping_timing(before, after, component):
    """Classify the first observed exposure; sparse histories can be ambiguous.

    Early support requires hydrogen exhaustion and central helium >= 0.9.
    This is a timing screen, not validation of an Ertl progenitor match.
    """
    prefix = f"S{component}_"
    hydrogen = number(after.get(prefix + "center_h1"))
    helium = number(after.get(prefix + "center_he4"))
    previous_helium = number(before.get(prefix + "center_he4"))
    if not (hydrogen <= 0.01 and math.isfinite(helium)):
        return "ambiguous"
    if helium >= 0.9:
        if (before.get(prefix + "state") == "H-rich_Shell_H_burning"
                and after.get("step_names") == "step_CE"):
            return "before_helium_ignition_supported"
        return "early_helium_supported"
    if helium <= 0.01 and previous_helium <= 0.01:
        return "after_helium_depletion_confirmed"
    if 0.01 < helium < 0.9 and 0.01 < previous_helium < 0.9:
        return "during_helium_burning_supported"
    return "ambiguous"


def summarize_history(history):
    """Retain resolved timing and mark unresolved intervals explicitly."""
    result = {}
    records = history.to_dict("records")
    interaction = next((r for r in records if
                        str(r.get("event", "")).startswith(("oCE", "oDoubleCE", "oRLO"))
                        or r.get("state") in ("RLO1", "RLO2", "contact")
                        or str(r.get("mass_transfer_case", "None"))
                        not in ("None", "nan", "<NA>", "", "no_RLO")), None)
    result["first_interaction_record_time_yr"] = (
        interaction.get("time") if interaction else None)
    result["first_interaction_record_event"] = (
        interaction.get("event") if interaction else None)
    for component in (1, 2):
        prefix = f"S{component}_"
        index = next((i for i, r in enumerate(records)
                      if str(r.get(prefix + "state", "")).startswith("stripped_He")), None)
        result[prefix + "stripping_timing"] = "not_observed"
        result[prefix + "first_strip_time_yr"] = None
        result[prefix + "first_strip_step"] = None
        result[prefix + "first_strip_mass"] = None
        result[prefix + "first_strip_center_he4"] = None
        if index is not None:
            after = records[index]
            before = records[max(0, index - 1)]
            result[prefix + "stripping_timing"] = stripping_timing(before, after, component)
            result[prefix + "first_strip_time_yr"] = after["time"]
            result[prefix + "first_strip_step"] = after.get("step_names")
            result[prefix + "first_strip_mass"] = after.get(prefix + "mass")
            result[prefix + "first_strip_center_he4"] = after.get(prefix + "center_he4")
    return result


def roche_fraction(mass, companion_mass):
    q_root = (mass / companion_mass) ** (1.0 / 3.0)
    return 0.49 * q_root**2 / (0.6 * q_root**2 + math.log1p(q_root))


def ce_energy_budget(mass1, mass2, core1, core2, radius1, radius2,
                     lambda1, lambda2, separation, alpha, double_ce,
                     gravitational_constant, solar_mass, solar_radius):
    """Evaluate the same signed-energy convention as POSYDON's CE step.

    Separations/radii are in solar radii; returned energies are in erg.
    A missing lambda is treated as zero binding energy by POSYDON, and flagged.
    """
    scale = gravitational_constant * solar_mass**2 / solar_radius
    missing_lambda = not math.isfinite(number(lambda1))
    binding = (0.0 if missing_lambda else
               scale * mass1 * (mass1 - core1) / (lambda1 * radius1))
    if double_ce:
        missing_lambda |= not math.isfinite(number(lambda2))
        if math.isfinite(number(lambda2)):
            binding += scale * mass2 * (mass2 - core2) / (lambda2 * radius2)
    initial_orbital_energy = -scale * mass1 * mass2 / (2 * separation)
    post_energy = initial_orbital_energy - binding / alpha
    predicted_separation = -scale * core1 * core2 / (2 * post_energy)
    return {
        "binding_energy_erg": binding,
        "initial_orbital_energy_erg": initial_orbital_energy,
        "fast_separation_rsun": predicted_separation,
        "missing_lambda_treated_as_zero_binding": bool(missing_lambda),
        "energy_scale_erg_rsun_per_msun2": scale,
    }


def fast_detachment_estimate(budget, core1, core2, radius1, radius2):
    """Necessary geometric energy test with fixed masses/radii.

    It ignores the subsequent stable-transfer/wind phase and radius responses
    to a different alpha. It does not predict the full CE-survival threshold.
    """
    if not all(math.isfinite(number(x)) and x > 0
               for x in (core1, core2, radius1, radius2)):
        return {}
    minimum_separation = max(radius1 / roche_fraction(core1, core2),
                             radius2 / roche_fraction(core2, core1))
    final_energy = (-budget["energy_scale_erg_rsun_per_msun2"]
                    * core1 * core2 / (2 * minimum_separation))
    available_energy = budget["initial_orbital_energy_erg"] - final_energy
    binding = budget["binding_energy_erg"]
    alpha_min = binding / available_energy if available_energy > 0 else math.nan
    return {
        "fast_min_detached_separation_rsun_estimate": minimum_separation,
        "fast_available_orbital_energy_erg_estimate": available_energy,
        "alpha_min_fast_detachment_estimate": alpha_min,
        "fast_separation_margin": budget["fast_separation_rsun"] / minimum_separation,
    }


def orbital_screen(binary, screen):
    """Screen the immediate post-collapse orbit, not present-day Gaia detection."""
    result = {"bound_ns_ms": False, "orbital_screen_pass": False,
              "periastron_detached": False, "M_NS": None, "M_companion": None}
    stars = (binary.star_1, binary.star_2)
    if stars[0].state == "NS" and stars[1].state == MS_STATE:
        ns, companion = stars
    elif stars[1].state == "NS" and stars[0].state == MS_STATE:
        companion, ns = stars
    else:
        return result
    result.update(M_NS=number(ns.mass), M_companion=number(companion.mass))
    eccentricity = number(binary.eccentricity)
    separation = number(binary.separation)
    if not (binary.state == "detached" and 0 <= eccentricity < 1 and separation > 0):
        return result
    result["bound_ns_ms"] = True
    companion_radius = 10 ** number(companion.log_R)
    companion_lobe = (separation * (1 - eccentricity)
                      * roche_fraction(companion.mass, ns.mass))
    result["periastron_detached"] = companion_radius < companion_lobe
    result["orbital_screen_pass"] = bool(
        result["periastron_detached"]
        and screen["companion_mass_min"] <= companion.mass <= screen["companion_mass_max"]
        and screen["period_min_days"] <= number(binary.orbital_period) <= screen["period_max_days"])
    return result


def make_diagnostic_ce_class():
    """Wrap the installed CE implementation; do not change its evolution."""
    from posydon.binary_evol.CE.step_CEE import StepCEE
    from posydon.utils import constants
    from posydon.utils.common_functions import orbital_separation_from_period

    class DiagnosticCEE(StepCEE):
        def __call__(self, binary):
            donor_number = 1 if str(binary.event).endswith("1") else 2
            self.record = {"donor": donor_number, "before": snapshot(binary)}
            self.budget = None
            self.fast_radii = None
            try:
                super().__call__(binary)
            except Exception as exc:
                self.record["error"] = f"{type(exc).__name__}: {exc}"
                raise
            finally:
                self.record["after"] = snapshot(binary)
                if self.budget is not None:
                    self.record.update(self.budget)
                    if self.fast_radii is not None:
                        self.record.update(fast_detachment_estimate(
                            self.budget, self.record["core1_mass"],
                            self.record["core2_mass"], *self.fast_radii))
                binary.study_ce_records.append(self.record.copy())

        def CEE_simple_alpha_prescription(self, binary, donor, comp_star,
                                         lambda1_CE, mc1_i, rc1_i, donor_type,
                                         lambda2_CE, mc2_i, rc2_i, comp_type,
                                         double_CE=False, **kwargs):
            alpha = self.common_envelope_efficiency
            self.record.update(alpha=alpha, lambda1=lambda1_CE, lambda2=lambda2_CE,
                               core1_mass=mc1_i, core2_mass=mc2_i,
                               core1_radius=rc1_i, core2_radius=rc2_i,
                               double_ce=double_CE,
                               core_h_fraction=self.core_definition_H_fraction,
                               post_ce_option=self.common_envelope_option_after_succ_CEE)
            separation = orbital_separation_from_period(
                binary.orbital_period, donor.mass, comp_star.mass)
            self.budget = ce_energy_budget(
                donor.mass, comp_star.mass, mc1_i, mc2_i,
                10**donor.log_R, 10**comp_star.log_R, lambda1_CE, lambda2_CE,
                separation, alpha, double_CE,
                constants.standard_cgrav, constants.Msun, constants.Rsun)
            self.fast_radii = (rc1_i, rc2_i)
            return super().CEE_simple_alpha_prescription(
                binary, donor, comp_star, lambda1_CE, mc1_i, rc1_i, donor_type,
                lambda2_CE, mc2_i, rc2_i, comp_type, double_CE=double_CE, **kwargs)

        def adjust_companion_radius(self, *args, **kwargs):
            radius = super().adjust_companion_radius(*args, **kwargs)
            if self.fast_radii is not None:
                self.fast_radii = (self.fast_radii[0], radius)
            return radius

        def CEE_adjust_post_CE_core_masses(self, *args, **kwargs):
            result = super().CEE_adjust_post_CE_core_masses(*args, **kwargs)
            self.fast_radii = (result[1], result[3])
            return result

    return DiagnosticCEE


def json_safe(value):
    """Write strict JSON without NaN or NumPy scalar values."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value
