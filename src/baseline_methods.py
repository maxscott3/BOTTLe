from typing import Callable, Tuple

import numpy as np


def euler_step(
        X: np.ndarray,  
        h: float, 
        stoich_matrix: np.ndarray,
        prop_funcs: Callable[[np.ndarray], np.ndarray],
        **kwargs
        ) -> np.ndarray:
    """Euler type tau-leaping step, see Definition 1"""
    prop_funcs_ = prop_funcs(X)
    poisson_rvs = np.random.poisson(lam= prop_funcs_ * h)
    X += poisson_rvs @ stoich_matrix
    X = np.clip(X, a_min = 0, a_max = None)
    return X, 0          #include 0 for degenerate event count


def thetatrap_step(
        X: np.ndarray, 
        h: float,
        stoich_matrix: np.ndarray, 
        prop_funcs: Callable[[np.ndarray], np.ndarray], 
        theta: float, 
        alpha1: float, 
        alpha2: float
        ) -> Tuple[np.ndarray, int]:
    """Theta trapezoidal method step, see Definition 2"""
    prop_funcs_1 = prop_funcs(X)
    poisson_rvs_1 = np.random.poisson(lam= prop_funcs_1 * h * theta)
    X_star = X + poisson_rvs_1 @ stoich_matrix
    X_star = np.clip(X_star, a_min = 0, a_max = None)
    prop_funcs_2 = prop_funcs(X_star)
    phi = (alpha1 * prop_funcs_2 - alpha2 * prop_funcs_1)
    idx_non_degenerate = np.where((phi >= 0))[0]
    poisson_rvs_2 = np.random.poisson(lam=phi[idx_non_degenerate] * (1 - theta) * h)
    X = X_star + poisson_rvs_2 @ stoich_matrix[idx_non_degenerate]
    degenerate_event_count = len(np.where(phi < 0)[0])
    X = np.clip(X, a_min = 0, a_max = None)
    return X, degenerate_event_count
