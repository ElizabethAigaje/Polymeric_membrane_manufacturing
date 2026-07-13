"""
Sensitivity analysis for membrane manufacturing LCA

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
import qsdsan as qs
import os

from _model_membrane_ import create_model


def sensitivity_analysis(sys, num_samples, membrane_area_per_year,
                         polymer_option, solvent_recycling,
                         threshold=0.1, save_path=None):
    """
    Perform Spearman correlation sensitivity analysis on membrane
    manufacturing LCA and plot a heatmap of correlations.

    Parameters
    ----------
    sys : qsdsan.System
        Simulated membrane system returned by system_membrane().
    num_samples : int
        Number of Latin Hypercube samples.
    membrane_area_per_year : float
        Annual membrane output [m²/year].
    polymer_option : str
        'PSF', 'CA', or 'CA_bioAA'.
    solvent_recycling : str
        'yes' or 'no'.
    threshold : float, optional
        Minimum absolute Spearman correlation to include a parameter
        in the heatmap. Default is 0.1 — parameters with all
        correlations below this are excluded for clarity.
    save_path : str, optional
        Folder to save the heatmap. Defaults to current script folder.

    Returns
    -------
    sensitivity_model : qsdsan.Model
        Model with all samples and results.
    r_df : pd.DataFrame
        Spearman correlation coefficients (parameters × indicators).
    p_df : pd.DataFrame
        P-values for each correlation.
    """

    save_path = save_path or os.path.dirname(__file__)

    # =========================================================================
    # 1. Create model with ONLY technological parameters
    # =========================================================================
    sensitivity_model = create_model(
        sys                    = sys,
        analysis               = 'technological',
        membrane_area_per_year = membrane_area_per_year,
        polymer_option         = polymer_option,
        solvent_recycling      = solvent_recycling,
    )

    np.random.seed(3221)

    samples = sensitivity_model.sample(N=num_samples, rule='L')
    sensitivity_model.load_samples(samples)
    sensitivity_model.evaluate()

    # Remove rows with NaN (failed simulations)
    sensitivity_model.table = sensitivity_model.table.dropna()
    n_valid = len(sensitivity_model.table)
    print(f'Valid samples: {n_valid} / {num_samples}')

    # =========================================================================
    # 2. Get Spearman correlations
    # =========================================================================
    r_df, p_df = qs.stats.get_correlations(sensitivity_model, kind='Spearman')

    # =========================================================================
    # 3. Process correlation dataframe for plotting
    # =========================================================================

    # --- Clean column names (metric names only, no units) ---
    # r_df columns are MultiIndex tuples: ('LCA', 'Global warming [kg CO2-Eq/m²]')
    # Extract just the metric name before the bracket
    clean_metric_names = [col[1].split(' [')[0] for col in r_df.columns]
    r_df.columns       = clean_metric_names

    # --- Clean row index (parameter names only, no units) ---
    # r_df index are MultiIndex tuples: ('DT101', 'Polymer fraction [w/w %/100]')
    clean_param_names = [idx[1].split(' [')[0] for idx in r_df.index]
    r_df.index        = clean_param_names

    # --- Filter parameters below threshold ---
    # Keep only parameters that have at least one correlation >= threshold
    mask        = r_df.abs().max(axis=1) >= threshold
    r_df_filtered = r_df[mask]

    n_removed = (~mask).sum()
    if n_removed > 0:
        print(f'\nRemoved {n_removed} parameters below threshold '
              f'(|rho| < {threshold}):')
        for p in r_df.index[~mask]:
            print(f'  {p}')

    print(f'\nParameters in heatmap: {len(r_df_filtered)}')
    print(f'Indicators in heatmap: {len(r_df_filtered.columns)}')


    # =========================================================================
    # 4. Print top correlations per indicator
    # =========================================================================
    print('\n=== Top 3 parameters per indicator ===')
    for indicator in r_df_filtered.columns:
        top3 = r_df_filtered[indicator].abs().nlargest(3)
        print(f'\n{indicator}:')
        for param, val in top3.items():
            sign = '+' if r_df_filtered.loc[param, indicator] > 0 else '-'
            print(f'  {sign}{val:.3f}  {param}')

    # =========================================================================
    # 5. Heatmap plot
    # =========================================================================
    n_params     = len(r_df_filtered)
    n_indicators = len(r_df_filtered.columns)

    plt.style.use('default')
    font = {'family': 'Arial', 'size': 8}
    plt.rc('font', **font)

    # Figure size scales with number of parameters
    cm_to_in   = 1 / 2.54
    fig_width  = max(12, n_indicators * 0.55) * cm_to_in
    fig_height = max(6, n_params * 0.65) * cm_to_in

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        data       = r_df_filtered.T,   # rows=indicators, cols=parameters
        vmin       = -1,
        vmax       = 1,
        center     = 0,
        cmap       = sns.diverging_palette(280, 220, s=90, l=70, as_cmap=True),  #h1 is for negatives and 220 is blue, #h2 is positive and 356 is pink, s is saturation (0=gray, 1--= very vivid color),  (0 dark and 100 white)
        linewidths = 0.3,
        linecolor  = 'white',
        ax         = ax,
        cbar_kws   = {'label': "Spearman's ρ", 'shrink': 0.6},
    )

    ax.set_xlabel('Parameter', fontsize=9)
    ax.set_ylabel('Impact indicator', fontsize=9)
    ax.set_title(
        f'Sensitivity analysis — {polymer_option}, '
        f'recycling={solvent_recycling}, n={n_valid}',
        fontsize=9, pad=10
    )

    # Rotate x-axis labels for readability
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center', fontsize=6)
    fig.subplots_adjust(bottom=0.35)   # add after tight_layout()
    ax.tick_params(axis='y', labelrotation=0,  labelsize=7)

    fig.tight_layout()
    plt.show()

    # =========================================================================
    # 6. Save outputs
    # =========================================================================
    os.makedirs(save_path, exist_ok=True)

    prefix = f'sensitivity_{polymer_option}_recycling{solvent_recycling}'

    # Save heatmap figure
    fig_path = os.path.join(save_path, f'{prefix}_heatmap.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f'\nSaved heatmap -> {fig_path}')

    # Save correlation table as CSV
    csv_path = os.path.join(save_path, f'{prefix}_correlations.csv')
    r_df_filtered.to_csv(csv_path)
    print(f'Saved correlations -> {csv_path}')

    # Save p-values table as CSV
    # p_df.columns = clean_metric_names[:len(p_df.columns)]
    # p_df.index   = clean_param_names[:len(p_df.index)]
    # pval_path    = os.path.join(save_path, f'{prefix}_pvalues.csv')
    # p_df.to_csv(pval_path)
    # print(f'Saved p-values -> {pval_path}')

    return sensitivity_model, r_df, p_df