import numpy as np

from utils import write_to_csv
    

def protein_aggregation(X_0: int, max_length: int, params: np.ndarray, file_name: str) -> None: 
    """Nucleation-elongation model of protein aggregation, i.e. proteins
        adding onto each other or breaking away. params is a 3d array for 
        nucleation/elongation/depolymerization. Produces reaction with (max_length) species 
        with (2*max_length - 2) reactions. initially fully concentrated (X_0) at length 1 polymer"""

    initial_species = np.zeros(max_length, dtype=np.int64)
    initial_species[0] = X_0

    stoich_matrix = np.zeros((2*max_length - 2, max_length), dtype=np.int8)
    stoich_matrix[0][0] = -2
    stoich_matrix[0][1] = 1

    for i in range(1, max_length-1):
        stoich_matrix[i][0] = -1
        stoich_matrix[i][i] = -1
        stoich_matrix[i][i+1] = 1
    
    for i in range(1, max_length):
        stoich_matrix[i+max_length-2][i] = -1
        stoich_matrix[i+max_length-2][i-1] = 1
        stoich_matrix[i+max_length-2][0] = 1
    
    param_list = np.empty((2*max_length - 2), dtype=np.float64)
    param_list[0] = params[0]
    param_list[1:max_length-1] = params[1]
    param_list[max_length-1:2*max_length-2] = params[2]

    write_to_csv(file_name, initial_species, param_list, stoich_matrix)



def morphogen_gradient(X_0: int, num_bins: int, params: np.ndarray, file_name: str) -> None:
    """Discretised 1d Morphogen Gradient model, proteins move up and down a chain, where there is
        initial production at an endpoint and the chance for degradation at any point on the 
        chain"""
    
    initial_species = np.zeros(num_bins, dtype=np.int64)
    initial_species[0] = X_0

    stoich_matrix = np.zeros((3*num_bins - 1, num_bins), dtype=np.int8)
    stoich_matrix[0][0] = 1

    for i in range(1, num_bins):
        stoich_matrix[i][i-1] = -1
        stoich_matrix[i][i] = 1

    for i in range(1, num_bins):
        stoich_matrix[i + num_bins - 1][i-1] = 1
        stoich_matrix[i + num_bins - 1][i] = -1
    
    for i in range(0, num_bins):
        stoich_matrix[i + 2*num_bins - 1][i] = -1

    param_list = np.empty((3*num_bins - 1))
    param_list[0] = params[0]
    param_list[1:num_bins] = params[1]
    param_list[num_bins:2*num_bins - 1] = params[1]
    param_list[2*num_bins-1:3*num_bins-1] = params[2]

    write_to_csv(file_name, initial_species, param_list, stoich_matrix)
