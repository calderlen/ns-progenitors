import p_nss


def gaia_dr3_detection_probability(ra, dec, period, eccentricity, parallax, M_NS, M_star, G_star, n_sample=1000,
):

    return p_nss.p_nss(ra, dec, period, eccentricity, parallax, M_NS, M_star, 999.0, G_star, n_sample=n_sample, mag="app")