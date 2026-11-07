import numpy as np


def sqrt_psd_matrix(matrix: np.ndarray) -> np.ndarray:
    """Return positive semi-definite matrix square root """
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    sqrt_eigenvalues = np.sqrt(np.clip(eigenvalues, a_min=0, a_max=None))
    sqrt_matrix = eigenvectors @ np.diag(sqrt_eigenvalues) @ (eigenvectors.T)
    return sqrt_matrix

