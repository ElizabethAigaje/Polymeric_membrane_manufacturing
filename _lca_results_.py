"""
LCA Results for membrane manufacturing system

LC systems analysis to identify green polymer materials and methods 
for membrane manufacturing.

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 04/27/2026
"""

# %%
# --- Imports ---
import os
import numpy as np
import pandas as pd
import qsdsan as qs

# %%
class LCAResults:
    """
    Computes and saves LCA results per m² membrane in three tables.

    Table 1 — Total impacts per m² (one row per indicator)
    Table 2 — Absolute contribution per source per m²
    Table 3 — Percentage contribution per source (% of total per indicator)

    Parameters
    ----------
    lca : qs.LCA
        The LCA object returned by system_membrane()
    membrane_area_per_year : float
        Total membrane area produced per year, [m²/year]. The LCA object returned by system_membrane()
    polymer_option : str
        'PSF', 'CA', or 'CA_bioAA' 
    solvent_recycling : str
        'yes' or 'no'
    save_path : str, optional
        Folder where CSV and Excel files will be saved.
        Defaults to the folder where this script lives.

    """

    def __init__(self, lca, membrane_area_per_year, polymer_option,
                 solvent_recycling, save_path=None):

        self.lca                    = lca
        self.membrane_area_per_year = membrane_area_per_year
        self.polymer_option         = polymer_option
        self.solvent_recycling      = solvent_recycling
        self.save_path              = save_path or os.path.dirname(__file__)

        # Build all three tables immediately on creation
        self.table1, self.table2, self.table3 = self._build_tables()

    # ------------------------------------------------------------------
    def _build_tables(self):
        lca  = self.lca
        area = self.membrane_area_per_year

        # --- Step 1: get raw annual impacts from qsdsan ---
        stream_raw   = lca.get_impact_table('Stream', annual=True)
        other_raw    = lca.get_impact_table('Other',  annual=True)
        total_annual = lca.get_total_impacts(annual=True)
        indicators   = list(total_annual.keys())

        # --- Step 2: Table 1 — total impacts per m² ---
        table1 = pd.DataFrame({
            'Total per year [unit/year]' : [total_annual[i] for i in indicators],
            'Total per m² [unit/m²]'    : [total_annual[i] / area for i in indicators],
            'Unit'                       : [
                qs.ImpactIndicator.get_indicator(i).unit
                if qs.ImpactIndicator.get_indicator(i) else ''   #In case an indicator is not registered, it puts an empety string instead of crushing (if..else)
                for i in indicators
            ],
        }, index=pd.Index(indicators, name='Indicator'))  #rows name

        # --- Step 3: helper to extract numeric columns from impact tables ---
        def extract_impact_cols(df, indicators):
            """
            Parse get_impact_table() output into a dict:
            {source_label: {indicator: annual_value}}
            Clean the rows and columns we do not need from get_impact_table
            """
            contributions = {}
            for source in df.index:
                if source == 'Sum':  #skip row 'Sum'
                    continue
                row = {}   # for each Row (source), it create an empty dictionary that loops over all 18 impacts
                for ind in indicators:
                    matching = [               #in the dataframe columns it will search for: start with indicator ID, has a '[' (units are inside) and the word ratio is not there (exclude ratio columns)
                        c for c in df.columns
                        if c.startswith(ind) and '[' in c and 'Ratio' not in c
                    ]
                    row[ind] = float(df.loc[source, matching[0]]) if matching else 0.0  #if it matches, add to the row (if nothing matches, adds 0.0)
                contributions[source] = row  #stores the complete row under the soruce name
            return contributions

        stream_contribs = extract_impact_cols(stream_raw, indicators)  #ditionaries of each stream and other with all the indicators. 'solvent': { 'GWP':  1200000.0,..}
        other_contribs  = extract_impact_cols(other_raw,  indicators)

        # --- Step 4: Names (labels) for each source ---
        stream_labels = {
            'solvent'               : 'NMP (solvent)',
            'additive'              : 'PEG (additive)',
            'water_boresol'         : 'Water in bore solution',
            'non_solvent'           : 'Water in coagulation bath',
            'rinsing_water'         : 'Water in rinsing',
            'glycerol_conditioning' : 'Glycerol (conditioning agent)',
            'water_conditioning'    : 'Water in conditioning',
            'diluted_wastewater'    : 'Wastewater to off-site WWT',
        }

        other_labels = {
            'PSF_item'              : f'Polymer ({self.polymer_option})',
            'CA_item'               : f'Polymer ({self.polymer_option})',
            'CA_bioAA_item'         : f'Polymer ({self.polymer_option})',
            'PSF_module_item'       : 'PSF in module housing',
            'epoxy_item'            : 'Epoxy resin in module housing',
            'water_wwt_item'        : 'Water in wastewater dilution',
            'NaOH_item'             : 'NaOH in regeneration',
            'ethanol_item'          : 'Ethanol in regeneration',
            'water_rsolution_item'  : 'Water in regeneration',
            'nitrogen_item'         : 'Nitrogen',
            'electricity_item'      : 'Electricity',
            'steam_item'            : 'Steam',
            'NMP_item'              : 'Makeup NMP (solvent)',
        }

        # --- Step 5: combine all sources, normalize per m² in a unified dict with new labels ---
        all_contribs = {}

        for raw_id, vals in stream_contribs.items():
            label = stream_labels.get(raw_id, raw_id)    #the default label is the second, inf not get the frist from stream_labels
            all_contribs[label] = {ind: v / area for ind, v in vals.items()}

        for raw_id, vals in other_contribs.items():
            label = other_labels.get(raw_id, raw_id)
            all_contribs[label] = {ind: v / area for ind, v in vals.items()}

        # --- Step 6: Table 2 — absolute contributions per m² --- (table.loc[Name of the row] select 1 row, table.columns are column headers, and table.index are the row labels, table.index.name is the name of the index column (1st), table.loc['TOtal, 'GWP] select a specific cell
        table2 = pd.DataFrame(all_contribs).T   #Dict into a dataframe, rows are soruces and columns are indciators
        table2.index.name = 'Source'
        table2 = table2.reindex(columns=indicators, fill_value=0.0)  #makes sure the columes are in the same order as indicator list and if any missing it will put 0.0
        table2.loc['TOTAL'] = table1['Total per m² [unit/m²]']  #add a total row at the end

        # --- Step 7: Table 3 — percentage contributions ---
        totals = table2.loc['TOTAL']
        table3 = table2.drop('TOTAL').div(totals) * 100   #drop row total, divide all the other rows existing in tabl2 to by totals and make it %
        table3.index.name = 'Source'
        table3.loc['TOTAL (%)'] = table3.sum()

        return table1, table2, table3

    # ------------------------------------------------------------------
    def show(self):
        """Print all three tables to the console."""
        print('\n' + '='*70)
        print(f'TABLE 1 — Total LCA impacts per m² '
              f'[{self.polymer_option}, recycling={self.solvent_recycling}]'), #name of the scenario
        print('='*70)
        print(self.table1.to_string())  #.to_string forec to print everyghitn regardless the console width (plain, no dataframe(df))

        print('\n' + '='*70)
        print('TABLE 2 — Absolute contribution per source [per m²]')
        print('='*70)
        print(self.table2.to_string())

        print('\n' + '='*70)
        print('TABLE 3 — Percentage contribution per source [%]')
        print('='*70)
        print(self.table3.to_string(float_format='{:.2f}'.format))

    # ------------------------------------------------------------------
    def save(self, filename_prefix=None):
        """
        Save all three tables to one Excel file (3 sheets) and
        three separate CSV files.

        Parameters
        ----------
        filename_prefix : str, optional
            Prefix for output filenames.
            Defaults to 'LCA_polymer_option_solvent_recycling'.
        """
        if filename_prefix is None:
            filename_prefix = (
                f'LCA_{self.polymer_option}_'
                f'recycling{self.solvent_recycling}'
            )

        os.makedirs(self.save_path, exist_ok=True)

        # Excel with 3 sheets
        xlsx_path = os.path.join(
            self.save_path, f'{filename_prefix}_results.xlsx'
        )
        with pd.ExcelWriter(xlsx_path) as writer:
            self.table1.to_excel(writer,
                sheet_name='Table1_Total_per_m2')
            self.table2.to_excel(writer,
                sheet_name='Table2_Contributions_per_m2')
            self.table3.to_excel(writer,
                sheet_name='Table3_Percent_contributions')
        print(f'Saved Excel  -> {xlsx_path}')

        # # Individual CSVs
        for i, table in enumerate([self.table1, self.table2, self.table3], 1):
            csv_path = os.path.join(
                self.save_path,
                f'{filename_prefix}_table{i}.csv'
            )
            table.to_csv(csv_path)
            print(f'Saved Table {i} CSV -> {csv_path}')