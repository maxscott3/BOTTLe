import numpy as np

from src.parameters import ModelParameters, ExperimentParameters, FigureParameters
from src.bottle_method import bottle_step
from src.baseline_methods import euler_step, thetatrap_step
from src.models import TauLeapingMethod, SSAMethod
from src.data import initialise_prop_functions

from data.utils import read_csv
from experiments.figures import FID_experiment_figure, FID_experiment_figure_w_degenerate_counts



X_0, param_list, stoich_matrix = read_csv("data/protein_aggregation.csv")
prop_funcs_SSA, prop_funcs = initialise_prop_functions(param_list, stoich_matrix)

modelParameters = ModelParameters(
    X_0 = X_0,
    stoich_matrix = stoich_matrix,
    prop_funcs = prop_funcs,
    high_dimensional = True
)

modelParametersSSA = ModelParameters(
    X_0 = X_0,
    stoich_matrix = stoich_matrix,
    prop_funcs = prop_funcs_SSA,
    high_dimensional = True
)

experimentParameters1 = ExperimentParameters(
    T = 1,
    step_sizes = 0.5 * np.array([0.001, 0.0008, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002]),
    N_samples = 1000,
    theta = 2/3
)

experimentParameters2 = ExperimentParameters(
    T = 1,
    step_sizes = np.array([0.001, 0.0008, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002]),
    N_samples = 1000,
    theta = 2/3
)

figureParameters = FigureParameters(
    filename = "experiments/experiment_2_FID_w_degenerate_count.png",
    n_methods = 3,
    method_labels = ["Tau-leaping", r"$\theta$-trapezoidal", r"B$\theta$TTLe"],
    method_colors = ["#D55E00", "#009E73", "#0072B2"],
    method_linestyles = ["dashed", "dashdot", "solid"],
    method_linewidths = [1.6, 1.8, 1.8]
)



def main():
    SSA_method = SSAMethod(modelParametersSSA, experimentParameters1)
    euler_method = TauLeapingMethod(euler_step, modelParameters, experimentParameters1)
    thetatrap_method = TauLeapingMethod(thetatrap_step, modelParameters, experimentParameters2)
    bottle_method = TauLeapingMethod(bottle_step, modelParameters, experimentParameters2)

    fid_metrics = np.empty(shape=(3, len(experimentParameters1.step_sizes)))
    degenerate_event_count = np.empty(shape=(3, len(experimentParameters1.step_sizes)))
    n_function_evaluations = np.empty(shape=(3, len(experimentParameters1.step_sizes)))

    print("Sampling SSA method")
    SSA_endpoints = SSA_method.multiple_sample()
    true_sample_mean = np.mean(SSA_endpoints, axis = 0)
    true_sample_cov = np.cov(SSA_endpoints, rowvar = False)

    print("Sampling Tau-leaping method")
    fid_metrics[0], degenerate_event_count[0] = euler_method.FID_experiment(true_sample_mean, true_sample_cov)
    n_function_evaluations[0] = experimentParameters1.T / experimentParameters1.step_sizes

    print("Sampling Theta-trapezoidal method")
    fid_metrics[1], degenerate_event_count[1] = thetatrap_method.FID_experiment(true_sample_mean, true_sample_cov)
    n_function_evaluations[1] = 2 * experimentParameters2.T / experimentParameters2.step_sizes

    print("Sampling Bottle method")
    fid_metrics[2], degenerate_event_count[2] = bottle_method.FID_experiment(true_sample_mean, true_sample_cov)
    n_function_evaluations[2] = 2 * experimentParameters2.T / experimentParameters2.step_sizes

    FID_experiment_figure_w_degenerate_counts(
        n_function_evaluations,
        fid_metrics,
        degenerate_event_count,
        figureParameters
    )
    print("Generated FID metric figure")


if __name__ == "__main__":
    main()
