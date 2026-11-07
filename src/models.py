from typing import Callable, Tuple
import numpy as np
from scipy.sparse import csr_matrix
from copy import deepcopy
from tqdm import tqdm, trange      

from src.parameters import ModelParameters, ExperimentParameters
from src.utils import sqrt_psd_matrix




class TauLeapingMethod:
    """Produce samples of a tau-leaping method given model and experiment parameters"""
    def __init__(self,
            method_step: Callable, 
            modelParameters: ModelParameters,
            experimentParameters: ExperimentParameters,
            ):
        self.method_step = method_step
        self.X_0 = modelParameters.X_0
        self.T = experimentParameters.T
        self.step_sizes = experimentParameters.step_sizes
        self.N_samples = experimentParameters.N_samples
        self.N_species = modelParameters.stoich_matrix.shape[1]

        self.theta = experimentParameters.theta
        self.alpha1 = 1 / (2 * self.theta * (1-self.theta))
        self.alpha2 = (self.theta ** 2 + (1 - self.theta) ** 2 ) * self.alpha1 

        #improves for high dimensional (num_species > 50) runtime 
        if modelParameters.high_dimensional:
            self.stoich_matrix = csr_matrix(modelParameters.stoich_matrix)
        else:
            self.stoich_matrix = modelParameters.stoich_matrix

        self.sample_params = {"stoich_matrix": self.stoich_matrix,
                                "prop_funcs": modelParameters.prop_funcs,
                                "theta": experimentParameters.theta,
                                "alpha1": self.alpha1,
                                "alpha2": self.alpha2}

        
    def single_sample(self, h: float) -> Tuple[np.ndarray, int]:
        """Generate a single sample run of experiment for a specific step size, returning 
            final species count and total number of degenerate events"""
        total_degenerate_event_count = 0
        steps = int(self.T / h)
        X = deepcopy(self.X_0)
        for _ in range(steps):
            X, degenerate_event_count = self.method_step(X, h, **self.sample_params)
            total_degenerate_event_count += degenerate_event_count
        return X, total_degenerate_event_count


    def multiple_sample(self, h: float) -> Tuple[np.ndarray, np.ndarray, int]:
        """Generate multiple samples of experiment for a specific step size, returning
            the mean and covariance of the final species counts along with  total number 
            of degenerate events"""
        total_degenerate_event_count = 0
        final_species_counts = np.empty((self.N_samples, self.N_species), dtype=np.int64)
        for i in range(self.N_samples):
            final_species_count, degenerate_event_count = self.single_sample(h)
            final_species_counts[i] = final_species_count
            total_degenerate_event_count += degenerate_event_count
        return final_species_counts, total_degenerate_event_count


    def FID_experiment(self, true_sample_mean: np.ndarray, true_sample_cov: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray]:
        """Measure Fréchet inception distance between method final species count and and assumed 
            exact methodology, here SSA. Returns metric for each step size in experiment along with
            average per sample number of degenerate events"""
        fid_metrics = np.empty(len(self.step_sizes), np.float64)
        avg_degenerate_event_count = np.empty(len(self.step_sizes), np.float64)
        true_sample_cov_sqrt = sqrt_psd_matrix(true_sample_cov)
        for i, h in enumerate(tqdm(self.step_sizes)):
            final_species_counts, total_degenerate_event_count = self.multiple_sample(h)        
            method_sample_mean = np.mean(final_species_counts, axis=0)
            method_sample_cov = np.cov(final_species_counts, rowvar=False)
            term = true_sample_cov_sqrt @ method_sample_cov @ true_sample_cov_sqrt
            sqrt_term = sqrt_psd_matrix(term)
            trace = np.trace(method_sample_cov + true_sample_cov - 2 * sqrt_term)
            mean_diff = true_sample_mean - method_sample_mean
            fid_metrics[i] = np.sqrt(np.dot(mean_diff, mean_diff) + trace)
            avg_degenerate_event_count[i] = total_degenerate_event_count / self.N_samples
        return fid_metrics, avg_degenerate_event_count




class SSAMethod:
    """Produce "exact" samples using the Stochastic Simulation Algorithm"""
    def __init__(self,
            modelParameters: ModelParameters,
            experimentParameters: ExperimentParameters
            ):
        self.X_0 = modelParameters.X_0
        self.stoich_matrix = modelParameters.stoich_matrix
        self.prop_funcs_SSA = modelParameters.prop_funcs
        self.T = experimentParameters.T
        self.N_samples = experimentParameters.N_samples
        self.N_reactions = self.stoich_matrix.shape[0]
        self.N_species = self.stoich_matrix.shape[1]
        self.max_count = 10e5              #hard-coded limit on steps in SSA sample


    def single_sample(self) -> np.ndarray:
        """Generate a single sample of SSA using the model and experiment parameters"""
        X = deepcopy(self.X_0)
        count = 0
        t = 0
        idx = np.arange(0, self.N_reactions)    #initially update all propensity functions
        prop_funcs_ = np.zeros(self.N_reactions, np.float64) 
        while (t < self.T) and (count < self.max_count):
            prop_funcs_ = self.prop_funcs_SSA(X, idx, prop_funcs_)
            t += np.log(1 / np.random.random()) / np.sum(prop_funcs_)
            probs = prop_funcs_ / np.sum(prop_funcs_)
            idx = np.array([np.random.choice(self.N_reactions, p=probs)])
            X += self.stoich_matrix[idx][0]
            count += 1
        if t < self.T: 
            #raise error if max_count reached
            raise ValueError(f"Max count reached on SSA sample. Final time t={t}.")
        return X

    def multiple_sample(self) -> np.ndarray:
        """Generate multiple samples of SSA using the model and experiment parameters"""
        final_species_counts = np.empty((self.N_samples, self.N_species), dtype=np.int32)
        for i in trange(self.N_samples):
            final_species_counts[i] = self.single_sample()
        return final_species_counts
    


