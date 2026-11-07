from typing import Callable, Tuple
import numpy as np


def bottle_step(
        X: np.ndarray, 
        h: float, 
        stoich_matrix: np.ndarray,
        prop_funcs: Callable[[np.ndarray], np.ndarray], 
        theta: float, 
        alpha1: float, 
        alpha2: float
    ) -> Tuple[np.ndarray, int, int]:
    """Bottle method step, see Definition 3"""
    prop_funcs_1 = prop_funcs(X)
    poisson_rvs_1 = np.random.poisson(lam= prop_funcs_1 * h * theta)
    X_star = X + poisson_rvs_1 @ stoich_matrix
    X_star = np.clip(X_star, a_min=0, a_max=None)   #truncate negative species counts to 0
    prop_funcs_2 = prop_funcs(X_star)
    phi = (alpha1 * prop_funcs_2 - alpha2 * prop_funcs_1)
    idx_non_degenerate = np.where((phi >= 0))[0]
    idx_degenerate = np.where((phi < 0))[0]
    probs = ( (1 - theta) / theta) * phi[idx_degenerate] / prop_funcs_1[idx_degenerate] + 1
    binomials = np.random.binomial(n = poisson_rvs_1[idx_degenerate], p = probs)
    poisson_rvs_2 = np.random.poisson(lam=phi[idx_non_degenerate] * (1 - theta) * h)
    X += (poisson_rvs_1[idx_non_degenerate] + poisson_rvs_2) @ stoich_matrix[idx_non_degenerate]
    X += binomials @ stoich_matrix[idx_degenerate]
    degenerate_event_count = len(np.where(phi < 0)[0])
    X = np.clip(X, a_min=0, a_max=None)
    return X, degenerate_event_count
