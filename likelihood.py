import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm


def object_loglikelihoods(observed, errors, model_masses, weights) -> np.ndarray:
    observed, errors, model_masses, weights = (
        np.asarray(values, dtype=float)
        for values in (observed, errors, model_masses, weights)
    )
    positive = weights > 0
    terms = norm.logpdf(
        observed[:, None],
        loc=model_masses[None, positive],
        scale=errors[:, None],
    )
    return logsumexp(terms + np.log(weights[positive])[None, :], axis=1)


def loglikelihood(observed, errors, model_masses, weights) -> float:
    return float(object_loglikelihoods(observed, errors, model_masses, weights).sum())
