from pathlib import Path

import numpy as np
import pandas as pd
import pyvo

TAP_URL = "https://dc.g-vo.org/tap"
TABLE = "gedr3mock.generated_data"
COLUMNS = ["source_id", "ra", "dec", "l", "b", "parallax", "age", "feh", "popid", "a_g", "phot_g_mean_mag", "visibility_periods_used", "random_index"]
MOCK_CSV = Path(__file__).resolve().parent / "outputs" / "gedr3mock_local.csv"


def query_local_mock(n_rows=500_000, max_distance_kpc=2.0, cache_path=MOCK_CSV, overwrite=False):
    cache_path = Path(cache_path)

    if cache_path.exists() and not overwrite:
        return pd.read_csv(cache_path)

    min_parallax_mas = 1.0 / max_distance_kpc
    columns = ", ".join(COLUMNS)

    query = f"""
        SELECT TOP {int(n_rows)} {columns}
        FROM {TABLE}
        WHERE parallax >= {min_parallax_mas} AND parallax > 0
        ORDER BY random_index
    """

    service = pyvo.dal.TAPService(TAP_URL)
    result = service.run_async(query, maxrec=int(n_rows))
    df = result.to_table().to_pandas()
    df["distance_kpc"] = 1.0 / df["parallax"]

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_path, index=False)

    return df


def sample_galactic_positions(df, n, seed=None, replace=None):
    if replace is None:
        replace = n > len(df)

    rng = np.random.default_rng(seed)
    indices = rng.choice(len(df), size=n, replace=replace)

    return df.iloc[indices].reset_index(drop=True)


if __name__ == "__main__":
    df = query_local_mock()
    print(df.head())
    print(f"rows: {len(df):,}")
