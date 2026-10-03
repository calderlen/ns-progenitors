from pathlib import Path

import pandas as pd


DATA_DIR = Path(__file__).resolve().parent / "data"


def load_gaia(path: str | Path = DATA_DIR / "gaia2024_analysis_master.csv") -> pd.DataFrame:
    columns = ["name", "gaia_dr3_id", "M2_Msun", "M2_err_Msun"]
    gaia = pd.read_csv(path, usecols=columns, dtype={"name": str, "gaia_dr3_id": str})
    return gaia[columns]


def load_theory(
    engine: str = "W18",
    mass_loss_sequence: str = "standard_Mdot",
    *,
    progenitor_path: str | Path = DATA_DIR / "ertl2020_progenitor_grid.csv",
    remnant_path: str | Path = DATA_DIR / "ertl2020_remnant_outputs_REQUIRED_TEMPLATE.csv",
) -> pd.DataFrame:
    progenitor_columns = [
        "model_id", "mass_loss_sequence", "MHe_initial_Msun", "MpreSN_Msun",
        "MZAMS_equiv_Msun", "MCO_Msun", "MFe_Msun",
    ]
    remnant_columns = ["model_id", "central_engine", "outcome_NS_BH", "M_NS_g_final_Msun"]
    progenitors = pd.read_csv(progenitor_path, usecols=progenitor_columns, dtype={"model_id": str})
    remnants = pd.read_csv(
        remnant_path, usecols=remnant_columns,
        dtype={"model_id": str, "central_engine": str, "outcome_NS_BH": str},
    )
    progenitors = progenitors.loc[progenitors["mass_loss_sequence"] == mass_loss_sequence]
    remnants = remnants.loc[remnants["central_engine"] == engine]
    theory = progenitors.merge(remnants, on="model_id", how="left")
    return theory.sort_values("MZAMS_equiv_Msun").reset_index(drop=True)
