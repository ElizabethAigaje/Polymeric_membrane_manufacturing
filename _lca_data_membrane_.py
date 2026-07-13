"""
LCA data for membrane manufacturing system

LC systems analysis to identify green polymer materials and methods
for membrane manufacturing.

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 04/23/2026

Uses brightway2.5 (bw2data / bw2io / bw2calc).
Uses ecoinvent 3.10 cutoff.

For PSF, CA, and CA-bioAA: CFs are injected directly from excel worksheet analysis.

For all other streams: CFs are calculated via bw2calc from ecoinvent 3.10.

Workflow
--------
FIRST RUN (or after any change to items/indicators):
    1. Import ecoinvent 3.10 once:
           uncomment import_ecoinvent() call at the bottom, run, re-comment.
    2. Run: save_cf_data()
       -> writes  data/cf_dct_membrane.pckl
                  data/cf_dct_membrane_summary.csv   (human-readable check)

EVERY SUBSEQUENT RUN (normal simulation):
    cf_dct = _load_lca_data_membrane()
"""

# %%
import os, sys, pickle
import numpy as np
import pandas as pd
import qsdsan as qs

# brightway2.5  (the three core packages — no "import brightway2")
import bw2data as bd
import bw2io   as bi
import bw2calc as bc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
c_path    = os.path.dirname(__file__)   #the path of the current script
data_path = os.path.join(c_path, 'data')
os.makedirs(data_path, exist_ok=True)

__all__ = (
    'import_ecoinvent',
    'create_indicators',
    'compute_cfs_from_ecoinvent',
    'get_cf_data',
    'save_cf_data',
    '_load_lca_data_membrane',
)

# %%
# ===========================================================================
# 1.  Project + database setup
# ===========================================================================
#
# Database names produced by bw2io.import_ecoinvent_release for v3.10 cutoff:
#   "ecoinvent-3.10-cutoff"    — technosphere (market activities etc.)
#   "ecoinvent-3.10-biosphere" — elementary flows
#
# LCIA methods are stored as tuples in bd.methods, e.g.:
#   ('ReCiPe 2016 v1.03, midpoint (H)', 'climate change', 'GWP1000')

#PROJECT_NAME = 'MembraneManufacturingLCA'
EI_DB_NAME   = 'ecoinvent-3.10-cutoff'
BS_DB_NAME   = 'ecoinvent-3.10-biosphere'

bd.projects.set_current('default')


def import_ecoinvent(username: str, password: str):
    """
    Import ecoinvent 3.10 cutoff into the brightway2.5 project.
    Run ONCE, then comment out.  Credentials are your ecoinvent account.
    """
    bi.import_ecoinvent_release(
        version      = '3.10',
        system_model = 'cutoff',
        username     = username,
        password     = password,
    )

if 'ecoinvent-3.10-cutoff' in bd.databases:    #Option A
    print('ecoinvent 3.10 is already present in the project')
else:
   import_ecoinvent('Shi_Research_Group', 'ShiGroup!')   #
        
eidb = bd.Database(EI_DB_NAME)
bsdb = bd.Database(BS_DB_NAME)



# %%
# ===========================================================================
# 2.  ReCiPe 2016 Midpoint (H) indicator definitions
#
#     Each entry: (qsdsan_ID, bw2_method_tuple, unit)
#
#     The bw2_method_tuple is the exact key in bd.methods for ecoinvent 3.10.
#     To inspect available methods after importing ecoinvent, run:
#       [m for m in bd.methods if 'ReCiPe 2016' in str(m) and 'midpoint' in str(m).lower()]
# ===========================================================================

INDICATORS = [
    # (qsdsan ID,   brightway2.5 method tuple ,                                                                                   unit)
    ('TAP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'acidification: terrestrial',              'terrestrial acidification potential (TAP)'),      'kg SO2-Eq'),
    ('GWP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'climate change',                         'global warming potential (GWP100)'),             'kg CO2-Eq'),
    ('FETP',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'ecotoxicity: freshwater',                'freshwater ecotoxicity potential (FETP)'),     'kg 1,4-DCB-Eq'),
    ('METP',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'ecotoxicity: marine',                    'marine ecotoxicity potential (METP)'),         'kg 1,4-DCB-Eq'),
    ('TETP',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'ecotoxicity: terrestrial',               'terrestrial ecotoxicity potential (TETP)'),    'kg 1,4-DCB-Eq'),
    ('FFP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'energy resources: non-renewable, fossil','fossil fuel potential (FFP)'),                    'kg oil-Eq'),
    ('FEP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'eutrophication: freshwater',             'freshwater eutrophication potential (FEP)'),      'kg P-Eq'),
    ('MEP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'eutrophication: marine',                 'marine eutrophication potential (MEP)'),          'kg N-Eq'),
    ('HTPc',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'human toxicity: carcinogenic',           'human toxicity potential (HTPc)'),               'kg 1,4-DCB-Eq'),
    ('HTPnc',   ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'human toxicity: non-carcinogenic',       'human toxicity potential (HTPnc)'),              'kg 1,4-DCB-Eq'),
    ('IRP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'ionising radiation',                     'ionising radiation potential (IRP)'),             'kBq Co-60-Eq'),
    ('LOP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'land use',                               'agricultural land occupation (LOP)'),            'm2a crop-Eq'),
    ('SOP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'material resources: metals/minerals',    'surplus ore potential (SOP)'),                   'kg Cu-Eq'),
    ('ODP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'ozone depletion',                        'ozone depletion potential (ODPinfinite)'),        'kg CFC-11-Eq'),
    ('PMFP',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'particulate matter formation',           'particulate matter formation potential (PMFP)'),  'kg PM2.5-Eq'),
    ('HOFP',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'photochemical oxidant formation: human health', 'photochemical oxidant formation potential: humans (HOFP)'), 'kg NOx-Eq'),
    ('EOFP',    ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'photochemical oxidant formation: terrestrial ecosystems', 'photochemical oxidant formation potential: ecosystems (EOFP)'), 'kg NOx-Eq'),
    ('WCP',     ('ecoinvent-3.10', 'ReCiPe 2016 v1.03, midpoint (H)',  'water use',                              'water consumption potential (WCP)'),              'm3'),
]

# Convenience lookups
IND_BW2_METHOD = {ind_id: method for ind_id, method, unit in INDICATORS}
IND_UNIT       = {ind_id: unit   for ind_id, method, unit in INDICATORS}

INDICATOR_IDS = [ind_id for ind_id, _, _ in INDICATORS]

def create_indicators(replace=True):
    """
    Register all 18 ReCiPe 2016 v1.03 Midpoint (H) indicators in qsdsan.
    Registers directly from the INDICATORS list above.

    Also checks that the bw2 method tuples exist in bd.methods and prints
    a warning for any that are missing (so you can fix the name).
    """
    for ind_id, bw2_method, unit in INDICATORS:
        existing = qs.ImpactIndicator.get_indicator(ind_id)
        if existing and replace:
            stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
            existing.deregister()
            sys.stdout = stdout
        if not qs.ImpactIndicator.get_indicator(ind_id):
            qs.ImpactIndicator(ID=ind_id, unit=unit)

    # Diagnostic: warn about any method tuple not found in the project
    missing = [
        (ind_id, bw2_method)
        for ind_id, bw2_method, unit in INDICATORS
        if bw2_method not in bd.methods
    ]
    if missing:
        print(
            '\nWARNING: The following LCIA method tuples were not found in '
            'bd.methods for this project.\nThey may have slightly different '
            'names in your ecoinvent 3.10 import.\n'
            'Fix: [m for m in bd.methods if "ReCiPe 2016" in str(m)]\n'
        )
        for ind_id, method in missing:
            print(f'  {ind_id}: {method}')

    return qs.ImpactIndicator.get_all_indicators()


# %%
# ===========================================================================
# 3.  Pre-computed CF scores (PSF, CA, bio-based acetic acid)
#     Taken directly from your literature table — per 1 kg of material.
#     Positive = burden (material consumed by the process).
# ===========================================================================

PRECOMPUTED_CFS = {

    # -------  Polysulfone (PSF)  -----------------------------------------------
    'PSF_item': {
        'TAP'     :  0.03291226,
        'GWP'     : 12.61042229,
        'FETP' :  0.46038306,
        'METP' :  0.60234816,
        'TETP' : 33.52388878,
        'FFP'     :  5.42904828,
        'FEP'     :  0.00278960,
        'MEP'     :  0.00032087,
        'HTPc'    : 1.834842940,
        'HTPnc'   : 11.18946937,
        'IRP'     :  0.48456118,
        'LOP'     :  0.17326557,
        'SOP'     :  0.29346329,
        'ODP'     :  2.5025e-06,
        'PMFP'    :  0.01446499,
        'HOFP'    :  0.02436003,
        'EOFP'    :  0.02689790,
        'WCP'     :  0.05716553,
    },

    # -------  Cellulose Acetate (CA)  ------------------------------------------
    'CA_item': {
        'TAP'     :  0.03899915,
        'GWP'     : 15.3744,
        'FETP' :  0.530925,   
        'METP' :  0.6723225,   
        'TETP' : 36.94345,
        'FFP'     :  5.41617,
        'FEP'     :  0.003506,
        'MEP'     :  0.000370413,
        'HTPc'    :  1.87275,
        'HTPnc'   : 12.80245,
        'IRP'     :  0.7244075,
        'LOP'     :  0.9139335,
        'SOP'     :  0.432526,
        'ODP'     :  2.94575e-06,
        'PMFP'    :  0.01929905,
        'HOFP'    :  0.03052855,
        'EOFP'    :  0.03394025,
        'WCP'     :  0.16403025,
    },

    # -------  CA with Bio-based Acetic Acid  -------------------------------------------
    'CA_bioAA_item': {
        'TAP'     :  0.03520208,
        'GWP'     :  13.69431775,
        'FETP' :  0.47711400,
        'METP' :  0.60208948,   
        'TETP' :  34.70232271,
        'FFP'     :  4.65928836,
        'FEP'     :  0.00314507,
        'MEP'     :  0.00033432,
        'HTPc'    :  1.58900571,
        'HTPnc'   :  10.84379758,
        'IRP'     :  0.71673273,
        'LOP'     :  0.89862827,
        'SOP'     :  0.41603289,
        'ODP'     :  0.00000285,
        'PMFP'    :  0.01722962,
        'HOFP'    :  0.02699908,
        'EOFP'    :  0.02993879,
        'WCP'     :  0.14368738,
    },
    # -------  Nitrogen — no representative ecoinvent activity  -------------
    'nitrogen_item': {ind_id: 0.0 for ind_id, _, _ in INDICATORS}
}

# Module housing is the same polymer — reference the same values
PRECOMPUTED_CFS['PSF_module_item'] = PRECOMPUTED_CFS['PSF_item']

# %%
# ===========================================================================
# 4.  Ecoinvent 3.10 activity lookup
# ===========================================================================

def _find_activity(name_fragment, location=None):
    """
    Searches in ecoinvent to find the exactly one activity that matches
    with what we are looking for
    Return one ecoinvent 3.10 activity whose name contains name_fragment
    (case-insensitive) and whose location matches exactly (if given).

    Raises LookupError with a diagnostic hint if 0 or >1 matches found.
    """
    results = [
        a for a in eidb
        if name_fragment.lower() in a['name'].lower()
        and (location is None or a.get('location', '') == location)   #providing or not the location, lower-- lowercase
    ]
    if len(results) == 0:                                             #if nothing is found
        raise LookupError(
            f'No activity found: name~"{name_fragment}", location="{location}".\n'
            f'Search with: [a["name"] for a in eidb '
            f'if "{name_fragment[:20].lower()}" in a["name"].lower()]'
        )
    if len(results) > 1:                                             #if multiple activities matched, need to be more specific with location for example
        names = '\n  '.join(
            f"{a['name']} | {a.get('location', '')}" for a in results
        )
        raise LookupError(
            f'Multiple matches for name~"{name_fragment}", location="{location}":\n  {names}'
        )
    return results[0]                     #Returns the exact activity needed


def _select_ei_activities():
    """
    Map each item_ID to its ecoinvent 3.10 Activity object.
    Edit name fragments / locations here if your database uses different names.

    Returns
    -------
    acts : dict  {item_ID: bw2data.Activity}
    """
    acts = {}

    # NMP — primary solvent
    acts['NMP_item']               = _find_activity('market for N-methyl-2-pyrrolidone',  'GLO')

    # PEG — pore-forming additive
    acts['PEG_item']               = _find_activity('market for triethylene glycol',      'RoW')  #proxy

    # Process water — multiple streams share same activity, kept separate for bookkeeping
    _water = _find_activity('market for water, deionised', 'RoW')
    acts['water_boresol_item']     = _water
    acts['nonsolvent_item']        = _water
    acts['water_rinsing_item']     = _water
    acts['water_rsolution_item']     = _water
    acts['water_conditioning_item']= _water
    
    acts['water_wwt_item']         = _find_activity('market for water, completely softened', 'US')

    # Glycerol — conditioning agent
    acts['glycerol_item']          = _find_activity('market for glycerine',                 'RoW')

    # NaOH — CA regeneration bath
    acts['NaOH_item']              = _find_activity('market for sodium hydroxide, without water', 'RoW')  #in 50% solution state', enough key to be unique

    # Ethanol — regeneration bath solvent (CA route only)
    acts['ethanol_item']           = _find_activity('market for ethanol, without water, in 95% solution state, from fermentation', 'RoW')

    # Epoxy resin — module assembly potting
    acts['epoxy_item']             = _find_activity('market for epoxy resin, liquid',       'RoW')

    #Outputs-waste
    # Wastewater treatment (sign is flipped in organize_cfs)
    acts['wastewater_item']        = _find_activity('treatment of wastewater, average, wastewater treatment', 'RoW')

    #Utilities
    # Electricity — 
    acts['electricity_item']       = _find_activity( 'market group for electricity, medium voltage', 'US')
    #Steam
    acts['steam_item']             = _find_activity('steam production, as energy carrier, in chemical industry', 'RoW')  #MJ; qsdsan carries in Kj -- untis conversion
    
    #In the LCA mian document: I should add to electricity (with the converison eqution) of heat exchangers
    # and the steam transform units from KJ (in qsdsan) to MJ (in the activity CF)

    return acts          

# %%
# ===========================================================================
# 5.  CF calculation via bw2calc
#     For each ecoinvent activity + each LCIA method: run a 1-unit LCA.
# ===========================================================================

def _compute_one_cf(activity, bw2_method_tuple):
    """
    Run the activity impact of 1 unit (FU) under bw2_method_tuple.
    Returns the LCIA score as a float, or np.nan if the method is missing.
    """
    if bw2_method_tuple not in bd.methods:  #So for example, _compute_one_cf(NMP_activity, GWP_tuple) returns something like 4.5 meaning 4.5 kg CO2-eq per kg of NMP.
        return np.nan
    lca = bc.LCA({activity: 1.0}, bw2_method_tuple)    #1 as it is per 1 FU and the tuple that has the indicator
    lca.lci()
    lca.lcia()
    return float(lca.score)


def compute_cfs_from_ecoinvent(ei_acts):
    """
    For every activity in ei_acts, compute CF scores for all INDICATORS.
    (Calls _compute_one_cf for every combination of material × indicator and 
    organises the answers into a dict.

    Returns
    -------
    cf_dct_ei : dict  {item_ID: {ind_ID: float}}
    """
    cf_dct_ei = {}
    n = len(ei_acts)
    for i, (item_id, activity) in enumerate(ei_acts.items(), 1):  
        print(
            f'  [{i}/{n}] {item_id}  ({activity["name"][:55]})',  #name of the activity and the ecoinvent backagroun: NMP_item and market for NMP ...
            flush=True,
        )
        scores = {
            ind_id: _compute_one_cf(activity, bw2_method)
            for ind_id, bw2_method, unit in INDICATORS
        }
        cf_dct_ei[item_id] = scores
    return cf_dct_ei


# %%
# ===========================================================================
# 6.  organize_cfs — sign conventions and unit conversions
# ===========================================================================

def organize_cfs(cf_dct_ei): ####something to fix here is the steam --- I have it as energy? no kg, what about utilities?
    """
    Combine ecoinvent-derived CFs with pre-computed literature CFs and
    apply sign / unit corrections.

    Sign convention:
      +  burden  (consumed input: NMP, PEG, electricity …)
      -  credit  (avoided burden or useful co-product)

    Wastewater:
      ecoinvent FU = 1 m3 treated (input to treatment).
      Our system produces wastewater as an output going to treatment.
      Therefore: flip sign AND convert m3 -> kg (divide by 997 kg/m3).

    Returns
    -------
    cf_dct : dict  {item_ID: {ind_ID: float}}
    """
    cf_dct = {}

    for item_id, scores in cf_dct_ei.items():
        if item_id == 'wastewater_item':
            cf_dct[item_id] = {k: -v / 997.0 for k, v in scores.items()}
        else:
            cf_dct[item_id] = dict(scores)   # positive burden, no conversion needed

    # Add pre-computed items (already per 1 kg, already correct sign)
    for item_id, scores in PRECOMPUTED_CFS.items():
        cf_dct[item_id] = dict(scores)

    return cf_dct   # all items including PSF, CA and CA(bio-aa) with the Ecoinvent based per FU 1 unit


# %%
# ===========================================================================
# 7.  create_items — register StreamImpactItems in qsdsan
# ===========================================================================

def create_items(cf_dct, replace=True):
    """
    Register one qs.StreamImpactItem per entry in cf_dct and populate CFs.
    This function takes your finished cf_dct (the dict with all CFs you just built) 
    and registers everything into qsdsan's memory so that when you later write 
    qs.SanStream('solvent', stream_impact_item=qs.StreamImpactItem.get_item('NMP_item')) 
    in your system file, qsdsan knows what NMP_item is and what its CFs are.
    Parameters
    ----------
    cf_dct  : dict  {item_ID: {ind_ID: float}}
    replace : bool  deregister and recreate if True (default)

    Returns
    -------
    items : list of qs.StreamImpactItem
    """
    items = []
    for item_id, scores in cf_dct.items():
        existing = qs.ImpactItem.get_item(item_id)   #get all the items in the directory we just built
        if existing and replace:
            existing.deregister()
        item = qs.StreamImpactItem(ID=item_id)   #creates a stream specific item
        for ind_id, cf_val in scores.items():
            if np.isnan(cf_val):   #previously if we had np so it does not crash
                continue
            if qs.ImpactIndicator.get_indicator(ind_id) is None:
                qs.ImpactIndicator(ID=ind_id, unit=IND_UNIT.get(ind_id, 'unknown'))
            item.add_indicator(ind_id, CF_value=cf_val)
        items.append(item)
    return items    #in qsdsan we will have the activities and the CF (populated)


# %%
# ===========================================================================
# 8.  get_cf_data — full pipeline (indicators -> activities -> CFs -> organize)
# ===========================================================================

def get_cf_data():
    """
    Run the complete CF generation pipeline. Returns cf_dct.
    This is the function called by save_cf_data().
    """
    print('Step 1/3  Registering impact indicators in qsdsan...')
    create_indicators()

    print('Step 2/3  Looking up ecoinvent 3.10 activities...')
    ei_acts = _select_ei_activities()

    print(
        f'Step 3/3  Computing CFs via bw2calc '
        f'({len(ei_acts)} activities x {len(INDICATORS)} indicators)...'
    )
    cf_dct_ei = compute_cfs_from_ecoinvent(ei_acts)

    print('         Organizing CFs and adding pre-computed items...')
    cf_dct = organize_cfs(cf_dct_ei)

    print('Done.')
    return cf_dct


# %%
# ===========================================================================
# 9.  save_cf_data — compute and persist to disk
# ===========================================================================

def save_cf_data():
    """
    Generate all CFs and save:
      data/cf_dct_membrane.pckl         <- loaded at runtime
      data/cf_dct_membrane_summary.csv  <- human-readable check
    """
    cf_dct = get_cf_data()   #run the previous function before qsdsan

    #Saves cf_dct as a binary pickle file. Think of pickle as taking a Python 
    # object (your dict with all the CFs) and freezing it to disk exactly as it is.
    #  The 'wb' means "write binary". This is the file that _load_lca_data_membrane() 
    # reads every time you run your simulation — so you only need to run the heavy bw2calc 
    # calculations once.
    
    pckl_path = os.path.join(data_path, 'cf_dct_membrane.pckl')
    with open(pckl_path, 'wb') as f:
        pickle.dump(cf_dct, f)
    
    #Creates a human-readable CSV so you can open it in Excel and visually check the CF values. 
    #The .T transposes the DataFrame so rows are items (NMP_item, PEG_item...) and columns are 
    #indicators (TAP, GWP...), which is easier to read. This file is just for your inspection 
    #The simulation never uses it.
    summary = pd.DataFrame(cf_dct).T
    summary.index.name = 'item_ID'
    csv_path = os.path.join(data_path, 'cf_dct_membrane_summary.csv')
    summary.to_csv(csv_path)

    print(f'\nSaved CF pickle  -> {pckl_path}')
    print(f'Saved CSV summary -> {csv_path}')
    return cf_dct


# %%
# ===========================================================================
# 10. _load_lca_data_membrane — called at runtime in your system script
# ===========================================================================

def _load_lca_data_membrane():
    """
    Load pre-computed CFs from disk and register all ImpactIndicators and
    StreamImpactItems in qsdsan.

    Call this once before creating the system, e.g.:

        from _lca_data_membrane import _load_lca_data_membrane
        cf_dct = _load_lca_data_membrane()

    Returns
    -------
    cf_dct : dict  {item_ID: {ind_ID: float}}
    """

    #Every time Python starts a new session, qsdsan's registry is completely empty — 
    # it lives only in memory, not on disk. So before loading any CFs, you need to re-register 
    # the 18 ReCiPe indicators (TAP, GWP, etc.) into qsdsan. The replace=False means "only register 
    # if not already there" — so if somehow they're already registered, don't touch them.
    create_indicators(replace=False)   # re-register indicators (in-memory only)

    #Checks that the pickle file actually exists on disk before trying to open it. If you forgot to 
    #run save_cf_data() first, or the path is wrong, it gives you a clear error message telling you 
    #exactly what to do instead of a cryptic Python crash.
    
    pckl_path = os.path.join(data_path, 'cf_dct_membrane.pckl')
    if not os.path.exists(pckl_path):
        raise FileNotFoundError(
            f'CF data not found at {pckl_path}.\n'
            'Generate it first by running:  save_cf_data()'
        )
    with open(pckl_path, 'rb') as f:  #opens it and read it in binary 
        cf_dct = pickle.load(f)     

    #Takes the loaded cf_dct and registers all the StreamImpactItem objects into qsdsan's registry
    #After this line, items like NMP_item, PSF_item, wastewater_item etc. are all live in memory 
    #and ready to be linked to streams in your system.
    
    create_items(cf_dct, replace=True)  #from the memory we created, we bring it and is ready to be used in qsdsan
    return cf_dct

#save_cf_data()    #uncomment just for the first time
# %%
# ===========================================================================
# 
#