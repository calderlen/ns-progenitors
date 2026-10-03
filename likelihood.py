import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm


def object_loglikelihoods(data, errors, model, weights) -> np.ndarray:
    # converting the inputs to numpy arrays
    data, errors, model, weights = (np.asarray(values, dtype=float) for values in (observed, errors, model_masses, weights))
    # selecting all positive weights
    positive = weights > 0
    # probability density of data D given model M 
    terms = norm.logpdf(data[:, None], loc=model[None, positive], scale=errors[:, None],)
    
    return logsumexp(terms + np.log(weights[positive])[None, :], axis=1)


def loglikelihood(data, errors, model, weights) -> float:
    return float(object_loglikelihoods(data, errors, model, weights).sum())
