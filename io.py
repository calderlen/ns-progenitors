from pathlib import Path

import pandas as pd

from ertl import load_theory


DATA_DIR = Path(__file__).resolve().parent / "data"


def load_gaia(path: str | Path = DATA_DIR / "gaia2024_analysis_master.csv") -> pd.DataFrame:
    columns = ["name", "gaia_dr3_id", "M2_Msun", "M2_err_Msun"]
    gaia = pd.read_csv(path, usecols=columns, dtype={"name": str, "gaia_dr3_id": str})
    return gaia[columns]
