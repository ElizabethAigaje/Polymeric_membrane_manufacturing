"""
Uncertainty analysis for membrane manufacturing LCA

LC systems analysis to identify green polymer materials and methods
for membrane manufacturing.

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 04/29/2026
"""

# --- Imports ---
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import qsdsan as qs
import os

from _model_membrane_ import create_model

def uncertainty_analysis(sys, analysis, num_samples, membrane_area_per_year,
                         polymer_option, solvent_recycling):
    """
    Perform uncertainty analysis on membrane manufacturing LCA and
    plot KDE curves for all 18 ReCiPe midpoint indicators.

    Parameters
    ----------
    sys : qsdsan.System
        Simulated membrane system returned by system_membrane().
    analysis : str
        Type of uncertainty to include: 'all', 'technological', or 'background'.
    num_samples : int
        Number of Latin Hypercube samples.
    membrane_area_per_year : float
        Annual membrane output [m²/year], used to normalize LCA metrics.
    polymer_option : str
        'PSF', 'CA', or 'CA_bioAA'.
    solvent_recycling : str
        'yes' or 'no'.

    Returns
    -------
    model_uncertainty : qsdsan.Model
        Model object with all samples and results stored in model_uncertainty.table.
    stats : dict
        Percentile statistics for each indicator.
    """

    # =========================================================================
    # 1. Create model and run sampling
    # =========================================================================
    model_uncertainty = create_model(
        sys                    = sys,
        analysis               = analysis,
        membrane_area_per_year = membrane_area_per_year,
        polymer_option         = polymer_option,
        solvent_recycling      = solvent_recycling,
    )

    np.random.seed(3221)

    samples = model_uncertainty.sample(N=num_samples, rule='L')
    model_uncertainty.load_samples(samples)
    model_uncertainty.evaluate()

    # =========================================================================
    # 2. Extract LCA results — FIX: filter to only actual metric columns
    # =========================================================================

    # Get the full LCA section of the table
    LCA_table = model_uncertainty.table.loc[:, 'LCA']

    # Build actual metric column names from model.metrics
    # This excludes parameter columns like 'NMP background [-]'
    # that qsdsan puts in the LCA section when element='LCA'
    actual_metric_names = [
        f'{m.name} [{m.units}]'
        for m in model_uncertainty.metrics
        if m.element == 'LCA'
    ]

    # Filter table to only the 18 actual metric columns
    LCA_results = LCA_table[actual_metric_names]
    LCA_metrics = actual_metric_names

    # =========================================================================
    # 3. Get baseline values — FIX: use tuple key access on MultiIndex
    # =========================================================================
    baseline_dict = {}
    for metric in LCA_metrics:
        try:
            baseline_dict[metric] = float(
                model_uncertainty.metrics_at_baseline()[('LCA', metric)]
            )
        except KeyError:
            # Fallback — use median if key not found
            baseline_dict[metric] = float(
                LCA_results[metric].replace([np.inf, -np.inf], np.nan).median()
            )
            print(f'WARNING: baseline not found for {metric}, using median')

    # Separate metric names and units for axis labels
    metric_names = [m.split(' [')[0]      for m in LCA_metrics]
    metric_units = [m.split(' [')[1][:-1] for m in LCA_metrics]

    # =========================================================================
    # 4. Calculate percentile statistics
    # =========================================================================
    percentiles = [2.5, 5, 25, 50, 75, 95, 97.5]
    stats = {}

    for metric in LCA_metrics:
        values = LCA_results[metric].replace([np.inf, -np.inf], np.nan).dropna()
        stats[metric] = dict(zip(
            [f'p{p}' for p in percentiles],
            np.percentile(values, percentiles)
        ))

    print('\n=== Uncertainty Analysis Results ===')
    for metric, pcts in stats.items():
        print(f'\n{metric}:')
        for p_label, v in pcts.items():
            print(f'  {p_label}: {v:.4e}')

    # =========================================================================
    # 5. KDE plots for all 18 indicators
    # =========================================================================
    n_indicators  = len(LCA_metrics)
    n_cols        = 3
    n_rows        = int(np.ceil(n_indicators / n_cols))

    plt.style.use('default')
    font = {'family': 'Arial', 'size': 8}
    plt.rc('font', **font)

    cm_to_in      = 1 / 2.54
    width_two_col = 17.1
    fig, axs = plt.subplots(n_rows, n_cols,
                             figsize=(width_two_col * cm_to_in,
                                      n_rows * 2.8 * cm_to_in))
    axs_flat = axs.flatten()

    colors = sns.color_palette(palette='husl', n_colors=n_indicators).as_hex()

    class MyScalarFormatter(ScalarFormatter):
        def _set_format(self):
            self.format = '%.2f'

    for i, (metric, color) in enumerate(zip(LCA_metrics, colors)):

        ax       = axs_flat[i]
        values   = LCA_results[metric].replace([np.inf, -np.inf], np.nan).dropna()
        baseline = baseline_dict[metric]   # ← uses dict, never baseline_values

        sns.kdeplot(values, ax=ax, fill=True, color=color,
                    bw_method='scott', linewidth=1.5)

        y_mid = 0.5 * (ax.get_ylim()[0] + ax.get_ylim()[1])
        ax.scatter(baseline, y_mid, marker='o', color=color,
                   s=15, zorder=5, label='Baseline')

        ax.axvline(stats[metric]['p5'],  color=color, linestyle='--',
                   linewidth=0.8, alpha=0.7, label='5th/95th pct')
        ax.axvline(stats[metric]['p95'], color=color, linestyle='--',
                   linewidth=0.8, alpha=0.7)

        ax.set_xlabel(f'{metric_names[i]}\n({metric_units[i]})',
                      ha='center', color=color, fontsize=7)
        ax.set_ylabel('')

        custom_formatter = MyScalarFormatter(useMathText=True)
        ax.xaxis.set_major_formatter(custom_formatter)
        ax.xaxis.major.formatter.set_powerlimits((0, 0))

        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.yaxis.set_ticks_position('none')
        ax.set_yticklabels([])
        ax.legend(fontsize=6)

    for j in range(n_indicators, len(axs_flat)):
        axs_flat[j].set_visible(False)

    fig.suptitle(
        f'LCA Uncertainty — {polymer_option}, recycling={solvent_recycling}, '
        f'n={num_samples}',
        fontsize=9, y=1.01
    )
    fig.tight_layout()
    plt.show()

    return model_uncertainty, stats