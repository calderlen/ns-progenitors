from pathlib import Path
import pandas as pd

THEORY_CSV = Path(__file__).resolve().parent / "data" / "ertl2020_remnants.csv"

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
        table = table.loc[
            table["mass_loss_sequence"] == mass_loss_sequence
        ]
    return table.reset_index(drop=True)
