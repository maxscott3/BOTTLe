from typing import Callable

import numpy as np
from copy import deepcopy



def initialise_prop_functions(param_list: np.ndarray, stoich_matrix: np.ndarray) \
        -> Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]:
    """Given propensity function parameters and stoichiometry matrix, efficiently define propensity
        functions. Avoids redundant updating of propensity functions with unchanged reactant counts"""
    
    if np.min(stoich_matrix) < -3:
        raise ValueError("Stoichiometry matrix invalid, highest order reactant too great")
    
    n_reactions = stoich_matrix.shape[0]
    linear, quad_1, quad_2, quad_double, cubic_1, cubic_2, cubic_3, cubic_sd_1, \
        cubic_sd_2, cubic_triple = [np.full(n_reactions, -1, dtype=np.int32) for _ in range(10)]
    reaction_products = np.full((n_reactions, n_reactions), False, dtype=np.int32)
    stoich_matrix_T = stoich_matrix.T
    reaction_types = np.zeros(n_reactions, dtype=np.int8)
    for i, stoich_vector in enumerate(stoich_matrix):
        n_reactants = np.sum(np.where(stoich_vector <= -1, -stoich_vector, 0))
        r_single= np.argwhere(stoich_vector == -1)
        r_double = np.argwhere(stoich_vector == -2)
        r_triple = np.argwhere(stoich_vector == -3)

        if n_reactants > 3:
            raise ValueError(f"Too many reactants in reaction {i}")
        
        elif n_reactants == 3:
            if len(r_triple) == 1:
                reaction_types[i] = 6
                cubic_triple[i] = r_triple[0][0]
            elif len(r_double) == 1:
                reaction_types[i] = 5
                cubic_sd_1[i] = r_single[0][0]
                cubic_sd_2[i] = r_double[0][0]
            else:
                reaction_types[i] = 4
                cubic_1[i] = r_single[0][0]
                cubic_2[i] = r_single[1][0]
                cubic_3[i] = r_single[2][0]

        elif n_reactants == 2:
            if len(r_double) == 1:
                reaction_types[i] = 3
                quad_double[i] = r_double[0][0]
            else:
                reaction_types[i] = 2
                quad_1[i] = r_single[0][0]
                quad_2[i] = r_single[1][0]

        elif n_reactants == 1:
            reaction_types[i] = 1
            linear[i] = r_single[0][0]
            
        #note that n_reactants == 0 accounted for already in reaction_types

        for j in np.argwhere(stoich_vector != 0):
            reaction_products[i] = reaction_products[i] | np.where(stoich_matrix_T[j] != 0, True, False)

    def prop_funcs_SSA(X: np.ndarray, idx_list: np.ndarray, prop_funcs_: np.ndarray) -> np.ndarray:
        """Propensity functions given indices of reaction rates to update. Note that this is done 
            inplace."""
        
        #consider both change in reactants and products
        all_idx = np.union1d(np.argwhere(np.sum(reaction_products[idx_list], axis=0) > 0), idx_list)
        types = reaction_types[all_idx]
        new_prop_funcs_ = deepcopy(prop_funcs_)

        m = types == 1
        if m.any():
            idx = all_idx[m]
            new_prop_funcs_[idx] = param_list[idx] * X[linear[idx]]

        m = types == 2
        if m.any():
            idx = all_idx[m]
            new_prop_funcs_[idx] = param_list[idx] * X[quad_1[idx]] * X[quad_2[idx]]

        m = types == 3
        if m.any():
            idx = all_idx[m]
            x_d = X[quad_double[idx]]
            new_prop_funcs_[idx] = param_list[idx] * x_d * (x_d - 1) / 2

        m = types == 4
        if m.any():
            idx = all_idx[m]
            new_prop_funcs_[idx] = param_list[idx] * X[cubic_1[idx]] * X[cubic_2[idx]] * X[cubic_3[idx]]

        m = types == 5
        if m.any():
            idx = all_idx[m]
            x_d = X[cubic_sd_2[idx]]
            new_prop_funcs_[idx] = param_list[idx] * X[cubic_sd_1[idx]] * x_d * (x_d - 1) / 2

        m = types == 6
        if m.any():
            idx = all_idx[m]
            x_t = X[cubic_triple[idx]]
            new_prop_funcs_[idx] = param_list[idx] * x_t * (x_t - 1) * (x_t - 2) / 6

        m = types == 0
        if m.any():
            idx = all_idx[m]
            new_prop_funcs_[idx] = param_list[idx]

        return new_prop_funcs_
    
    def prop_funcs(X: np.ndarray):
        return prop_funcs_SSA(X, np.arange(0, n_reactions), np.zeros(n_reactions))
    
    return prop_funcs_SSA, prop_funcs
