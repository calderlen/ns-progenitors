from contextlib import chdir
from pathlib import Path
import sys


SELECTION_DIR = Path(__file__).resolve().parent / "dr3_nss_ast_orbits_selection"
sys.path.insert(0, str(SELECTION_DIR / "code"))

from p_nss import p_nss


def gaia_dr3_detection_probability(
    ra, dec, period, eccentricity, parallax, M_NS, M_star, G_star, n_sample=1000,
):

    # The upstream model reads its calibration files from the working directory.
    with chdir(SELECTION_DIR):
        return p_nss(
            ra, dec, period, eccentricity, parallax,
            M_NS, M_star, 999.0, G_star, n_sample=n_sample, mag="app",
        )
