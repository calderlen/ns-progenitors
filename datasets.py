from pathlib import Path

import pandas as pd


DATA_DIR = Path(__file__).resolve().parent / "data"
THEORY_CSV = DATA_DIR / "ertl2020" / "ertl2020_remnants.csv"
GAIA_CSV = DATA_DIR / "gaia2024" / "gaia2024_analysis_master.csv"


def load_theory(
    engine: str | None = None,
    mass_loss_sequence: str | None = None,
    *,
    path: str | Path = THEORY_CSV,
) -> pd.DataFrame:
    table = pd.read_csv(path)
    if engine is not None:
        table = table.loc[table["central_engine"] == engine]
    if mass_loss_sequence is not None:
        table = table.loc[table["mass_loss_sequence"] == mass_loss_sequence]
    return table.reset_index(drop=True)


def load_gaia(path: str | Path = GAIA_CSV) -> pd.DataFrame:
    columns = ["name", "gaia_dr3_id", "M2_Msun", "M2_err_Msun"]
    gaia = pd.read_csv(path, usecols=columns, dtype={"name": str, "gaia_dr3_id": str})
    return gaia[columns]
