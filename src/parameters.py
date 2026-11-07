from typing import Callable, List
import numpy as np
from dataclasses import dataclass


@dataclass
class ModelParameters:
    X_0: np.ndarray
    stoich_matrix: np.ndarray
    prop_funcs: Callable
    high_dimensional: bool


@dataclass
class ExperimentParameters:
    T: float
    step_sizes: np.ndarray
    N_samples: int
    theta: float


@dataclass
class FigureParameters:
    filename: str
    n_methods: int
    method_labels: List[str]
    method_colors: List[str]
    method_linestyles: List[str]
    method_linewidths: List[float]
    