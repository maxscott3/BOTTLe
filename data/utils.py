from typing import Tuple
import numpy as np
import pandas as pd


def write_to_csv(file_name: str, X_0: np.ndarray, param_list: np.ndarray, stoich_matrix: np.ndarray) \
        -> None:
    """Write initial species population, propensity function parameters and stoichiometry matrix 
        to a .csv file"""
    if file_name[-4:] != ".csv":
        file_name += ".csv"
    array = np.empty(shape=(len(param_list)+1, len(X_0)+1))
    array[0, 1:] = X_0
    array[1:, 0] = param_list
    array[1:, 1:] = stoich_matrix
    df = pd.DataFrame(array)
    df.to_csv(file_name, header=False, index=False)


def read_csv(file_name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract initial species population, propensity function parameters and stoichiometry matrix 
        from a .csv file"""
    model_data = pd.read_csv(file_name, header=None)
    X_0 = model_data.iloc[0, 1:].to_numpy(dtype=np.int64)
    param_list = model_data.iloc[1:, 0].to_numpy(dtype=np.float32)
    stoich_matrix = model_data.iloc[1:, 1:].to_numpy(dtype=np.int8)
    return X_0, param_list, stoich_matrix
