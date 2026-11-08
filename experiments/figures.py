import numpy as np
import matplotlib.pyplot as plt

from src.parameters import FigureParameters

def FID_experiment_figure(
        n_function_evaluations: np.ndarray,
        fid_metrics: np.ndarray,
        figureParameters: FigureParameters
        ):
    """Produce figure for FID experiment"""
    plt.figure(figsize=(8,5))
    for i in range(figureParameters.n_methods):
        plt.loglog(n_function_evaluations[i], 
                 fid_metrics[i], 
                 color = figureParameters.method_colors[i],
                 linestyle = figureParameters.method_linestyles[i],
                 linewidth = figureParameters.method_linewidths[i],
                 label = figureParameters.method_labels[i])
            
    plt.xlabel('number of propensity function evaluations', fontsize=12)
    plt.ylabel('FID metric', fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig(figureParameters.filename, bbox_inches='tight')


def FID_experiment_figure_w_degenerate_counts(
        n_function_evaluations: np.ndarray,
        fid_metrics: np.ndarray,
        degenerate_event_counts: np.ndarray,
        figureParameters: FigureParameters
        ):
    """Produce figure for FID experiment"""
    plt.figure(figsize=(8,8))
    fig, (ax0, ax1) = plt.subplots(
        2, 1, 
        height_ratios=[10, 7], 
        sharex=True,
        gridspec_kw={"hspace": 0.05})
    for i in range(figureParameters.n_methods):
        ax0.loglog(n_function_evaluations[i], 
                 fid_metrics[i], 
                 color = figureParameters.method_colors[i],
                 linestyle = figureParameters.method_linestyles[i],
                 linewidth = figureParameters.method_linewidths[i],
                 label = figureParameters.method_labels[i])
        if np.sum(degenerate_event_counts) > 0:
            ax1.loglog(n_function_evaluations[i], 
                 degenerate_event_counts[i], 
                 color = figureParameters.method_colors[i],
                 linestyle = figureParameters.method_linestyles[i],
                 linewidth = figureParameters.method_linewidths[i],
                 label = figureParameters.method_labels[i])
            
    ax0.set_xlabel('number of propensity function evaluations', fontsize=12)
    ax1.set_xlabel('number of propensity function evaluations', fontsize=12)
    ax0.set_ylabel('FID metric', fontsize=12)
    ax1.set_ylabel('degenerate event count', fontsize=12)
    fig.align_ylabels()
    ax0.tick_params(labelbottom=False)
    ax0.legend(fontsize=12)
    plt.savefig(figureParameters.filename, bbox_inches='tight')


