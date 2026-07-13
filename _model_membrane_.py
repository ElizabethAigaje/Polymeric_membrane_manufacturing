"""
Uncertainty and sensitivity analysis model for membrane manufacturing LCA

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
from chaospy import distributions as shape
import qsdsan as qs
import numpy as np
import pandas as pd
import os

from _lca_data_membrane_ import _load_lca_data_membrane, INDICATORS, INDICATOR_IDS

__all__ = ('create_model',)


def create_model(sys, analysis,membrane_area_per_year, polymer_option, solvent_recycling):
    """
    Create a qsdsan Model for uncertainty and sensitivity analysis.

    Parameters
    ----------
    sys : qsdsan.System
        The simulated membrane system returned by system_membrane().
    membrane_area_per_year : float
        Annual membrane output [m²/year], used to normalize LCA metrics.
    polymer_option : str
        'PSF', 'CA', or 'CA_bioAA', help to determine parameters in the model
    solvent_recycling : str
        'yes' or 'no'
    
    Returns
    -------
    model : qsdsan.Model
        Model object with parameters and metrics defined, ready for sampling.
    """

    # --- Initialize model ---
    model = qs.Model(sys)
    param  = model.parameter
    metric = model.metric

    # --- Get LCA object and flowsheet shortcuts ---
    lca      = sys.LCA
    flowsheet = qs.Flowsheet.flowsheet.default
    fs_stream = flowsheet.stream
    fs_unit   = flowsheet.unit

    # =========================================================================
    # PARAMETERS
    # =========================================================================
    # Parameters tell the model WHAT to vary and by HOW MUCH.
    # Each parameter needs:
    #   name        : readable label for plots and tables
    #   element     : which part of the system it belongs to (unit ID or 'LCA') This is just a label for organzing results
    #   kind        : 'coupled' (affects simulation) or 'isolated' (post-sim adjustment)
    #   units       : string for axis labels
    #   baseline    : the default/central value
    #   distribution: probability distribution describing uncertainty range
    #    -------------------------------------------------------------------------

    # Abreviations: ps = polymer solution, bs = bore solution, rs = regenerant solution, cs = conditioning solution, ww = wastewater

    def set_technological_params():
        """
        Parameters that directly affect the simulation outputs.
       
        """
        # --- System-level variable (aggregated and updated in MmebraneParams) --- (7)

        params= fs_unit.DT101.params   #It can be through any unit that has the params term, params is the same object across all units

        baseline = params.polymer_fraction  
        dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Polymer fraction', element='DT101', kind='coupled', units='w/w %/100', baseline=baseline, distribution=dist)
        def set_polymer_fraction(i):
            params.polymer_fraction = i
            params.refresh()  #call every time a parameter change during uncertainty to recalculate
        
        baseline = params.additive_fraction  
        dist = shape.Triangle(lower=0.04, midpoint=baseline, upper=0.16)
        @param(name='Additive fraction', element='M102', kind='coupled', units='w/w %/100', baseline=baseline, distribution=dist)
        def set_additive_fraction(i):
            params.additive_fraction = i
            params.refresh()
       
        dist=shape.Uniform(lower=0.8, upper=1.29)
        @param(name='Volume ratio: bs to ps', element='M101', kind ='coupled', units='-', baseline=params.ratio_boresol_to_polymersol, distribution=dist)  
        def set_ratio_boresol_to_polymersol(i):  
            params.ratio_boresol_to_polymersol=i 
            params.refresh()
        
        baseline = params.solventbore_fraction 
        dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Solvent fraction, bs.', element='M101', kind='coupled', units='w/w %/100', baseline=baseline, distribution=dist)
        def set_solventbore_fraction (i):
            params.solventbore_fraction  = i
            params.refresh()

        dist=shape.Uniform(lower=0.15, upper=0.35)
        @param(name='Glycerol concentration', element='M301', kind ='coupled', units='w/w %/100', baseline=params.conditioning_concentration, distribution=dist)  
        def set_conditioning_concentration(i):  
            params.conditioning_concentration=i 
            params.refresh()

        #For the heat exchangers, the HX units do not have the params attribute, 
        # so once the value is update through params (and with that all the other parameters that depend on T inside params); 
        # however the Hx.T is not updated, so we need to manually tell HX to have a new temperature 

        baseline = params.T_polymersol  
        dist = shape.Triangle(lower=40, midpoint=baseline, upper=60)
        @param(name='Dissolution temperature', element='HX101', kind='coupled', units='°C', baseline=baseline, distribution=dist)
        def set_T_polymersol(i):
            params.T_polymersol = i   #params
            params.refresh()
            # Also update HX101 and HX102 directly since they stored the value as a float
            # Note: HXutility expects temperature in Kelvin
            fs_unit.HX101.T = i + 273.15   # both temperatures are updated similarly
            fs_unit.HX102.T = i + 273.15   

        baseline = params.T_boresol  
        dist = shape.Triangle(lower=40, midpoint=baseline, upper=60)
        @param(name='Bore solution temperature', element='HX103', kind='coupled', units='°C', baseline=baseline, distribution=dist)
        def set_T_boresol(i):
            params.T_boresol = i
            params.refresh()

            fs_unit.HX103.T = i + 273.15   
        
        # --- S1: Bore and polymer solution preparation --- (10)
        
        baseline = fs_unit.DT101.tau
        dist = shape.Triangle(lower=6, midpoint=baseline, upper=35)
        @param(name='Dissolution time', element='DT101', kind='coupled', units='h', baseline=baseline, distribution=dist)
        def set_Dissolution_time(i):
            fs_unit.DT101.tau = i
        
        baseline = fs_unit.DT101.kW_per_m3
        dist = shape.Triangle(lower=12.6, midpoint=baseline, upper=22)
        @param(name='Dissolution energy', element='DT101', kind='coupled', units='kW/m^3', baseline=baseline, distribution=dist)
        def set_Dissolution_energy(i):
            fs_unit.DT101.kW_per_m3 = i
        
        baseline = fs_unit.DG101.tau
        dist = shape.Triangle(lower=2, midpoint=baseline, upper=24)
        @param(name='Degassing time, ps', element='DG101', kind='coupled', units='h', baseline=baseline, distribution=dist)
        def set_Degassing_tau_ps(i):
            fs_unit.DG101.tau = i
        
        dist = shape.Uniform(lower=0.04, upper=0.1)
        @param(name='Degassing energy, ps.', element='DG101', kind='coupled', units='kW/m^3', baseline=fs_unit.DG101.kW_per_m3, distribution=dist)
        def set_Degassing_energy_ps(i):
            fs_unit.DG101.kW_per_m3 = i
        
        dist = shape.Uniform(lower=700, upper=760)
        @param(name='Degassing vacuum pressure, ps', element='DG101', kind='coupled', units='torr', baseline=fs_unit.DG101.P, distribution=dist)
        def set_Degassing_vacumm_pressure(i):
            fs_unit.DG101.P = i
        
        dist = shape.Uniform(lower=0.9, upper=1.1)
        @param(name='Mixing time, bs', element='DT102', kind='coupled', units='h', baseline=fs_unit.DT102.tau, distribution=dist)
        def set_Mixing_time_bs(i):
            fs_unit.DT102.tau = i
        
        dist = shape.Uniform(lower=1, upper=1.5)
        @param(name='Mixing energy, bs', element='DT102', kind='coupled', units='kW/m^3', baseline=fs_unit.DT102.kW_per_m3, distribution=dist)
        def set_Mixing_energy_bs(i):
            fs_unit.DT102.kW_per_m3 = i
        
        baseline = fs_unit.DG102.nitrogen_demand
        dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Nitrogen demand', element='DG102', kind='coupled', units='L/m^2 membrane', baseline=baseline, distribution=dist)
        def set_Nitrogen_demand(i):
            fs_unit.DG102.nitrogen_demand = i
        
        dist = shape.Uniform(lower=0.9, upper=1.1)
        @param(name='Degassing time, bs', element='DG102', kind='coupled', units='h', baseline=fs_unit.DG102.tau, distribution=dist)
        def set_Degassing_tau_bs(i):
            fs_unit.DG102.tau = i
        
        dist = shape.Uniform(lower=0.04, upper=0.1)
        @param(name='Degassing energy, bs', element='DG102', kind='coupled', units='kW/m^3', baseline=fs_unit.DG102.kW_per_m3, distribution=dist)
        def set_Degassing_energy_bs(i):
            fs_unit.DG102.kW_per_m3 = i
        
        # --- S2: Extrusion, coagulation and rinsing --- (12)

        baseline = fs_unit.EX201.power_demand
        dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Extrusion power demand', element='EX201', kind='coupled', units='kW', baseline=baseline, distribution=dist)
        def set_Extrusion_power(i):
            fs_unit.EX201.power_demand = i
        
        if solvent_recycling== 'no':  #In solvent recycling we do not trest this variable, as it is valid only under especifc water use conditions (scenario)
            dist = shape.Uniform(lower=50, upper=250)
            @param(name='Non solvent demand', element='M201', kind='coupled', units='kg/m^2 membrane', baseline=fs_unit.M201.nonsolvent_per_membranearea, distribution=dist)
            def set_Nonsolvent_demand(i):
                fs_unit.M201.nonsolvent_per_membranearea = i
        
        baseline = fs_unit.HX201.T
        dist = shape.Triangle(lower=283.15, midpoint=baseline, upper=303.15)
        @param(name='Coagulation temperature', element='HX201', kind='coupled', units='K', baseline=baseline, distribution=dist)
        def set_Coagulation_temperature(i):
            fs_unit.HX201.T = i
        
        if polymer_option in ('CA', 'CA_bioAA'):
            
            baseline = fs_unit.RGT201.naoh_concentration
            dist = shape.Triangle(lower=0.05, midpoint=baseline, upper=0.2)
            @param(name='NaOH concentration', element='RGT201', kind='coupled', units='M', baseline=baseline, distribution=dist)
            def set_NaOH_concentration(i):
                fs_unit.RGT201.naoh_concentration = i
            
            baseline = fs_unit.RGT201.tau
            dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
            @param(name='Mixing time, rs', element='RGT201', kind='coupled', units='h', baseline=baseline, distribution=dist)
            def set_Mixing_time_rs(i):
                fs_unit.RGT201.tau = i
        
            baseline = fs_unit.RGT201.kW_per_m3
            dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
            @param(name='Mixing energy, rs', element='RGT201', kind='coupled', units='kW/m^3', baseline=baseline, distribution=dist)
            def set_Mixing_energy_rs(i):
                fs_unit.RGT201.kW_per_m3 = i
            
            baseline = fs_unit.RB201.tau
            dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
            @param(name='Regeneration time', element='RB201', kind='coupled', units='h', baseline=baseline, distribution=dist)
            def set_Regeneration_time(i):
                fs_unit.RB201.tau = i
            
            dist = shape.Uniform(lower=1, upper=1.5)
            @param(name='Regeneration energy demand', element='RB201', kind='coupled', units='kW/m^3', baseline=fs_unit.RB201.kW_per_m3, distribution=dist)
            def set_Regeneration_energy_demand(i):
                fs_unit.RB201.kW_per_m3 = i
            
            # dist = shape.Uniform(lower=0.01, upper=0.05)
            # @param(name='Purge fraction', element='SP201', kind='coupled', units='%/100', baseline=float(fs_unit.SP201.split), distribution=dist)
            # def set_Purge_fraction(i):
            #     fs_unit.SP201.split[:] = i

        if solvent_recycling== 'no':
            dist = shape.Uniform(lower=50, upper=250)
            @param(name='Rinsing water demand', element='M202', kind='coupled', units='kg/m^2 membrane', baseline=fs_unit.M202.rinsing_per_membranearea, distribution=dist)
            def set_Rinsing_water(i):
                fs_unit.M202.rinsing_per_membranearea = i

        dist = shape.Uniform(lower=298.15, upper=323.15)
        @param(name='Rinsing temperature', element='HX202', kind='coupled', units='K', baseline=fs_unit.HX202.T, distribution=dist)
        def set_Rinsing_temperature(i):
            fs_unit.HX202.T = i

        dist = shape.Uniform(lower=0.08, upper=1)
        @param(name='Rinsing time', element='RT201', kind='coupled', units='h', baseline=fs_unit.RT201.tau, distribution=dist)
        def set_Rinsing_time(i):
            fs_unit.RT201.tau = i

        # --- S3: Conditioning and drying --- (10)

        dist = shape.Uniform(lower=1, upper=1.21)
        @param(name='Glycerol pore filling fraction ', element='M301', kind='coupled', units='%/100', baseline=fs_unit.M301.fraction_filled, distribution=dist)
        def set_Glycerol_porefraction(i):
            fs_unit.M301.fraction_filled = i

        dist = shape.Uniform(lower=1, upper=1.21)
        @param(name='Water pore filling fraction', element='M302', kind='coupled', units='%/100', baseline=fs_unit.M302.fraction_filled, distribution=dist)
        def set_Water_porefraction(i):
            fs_unit.M302.fraction_filled = i

        baseline = fs_unit.MT301.tau
        dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Mixing time, cs', element='MT301', kind='coupled', units='h', baseline=baseline, distribution=dist)
        def set_Mixing_time_cs(i):
            fs_unit.MT301.tau = i
        
        dist = shape.Uniform(lower=1, upper=1.5)
        @param(name='Mixing energy, cs', element='MT301', kind='coupled', units='kW/m^3', baseline=fs_unit.MT301.kW_per_m3, distribution=dist)
        def set_Mixing_energy_cs(i):
            fs_unit.MT301.kW_per_m3 = i

        baseline = fs_unit.HX301.T
        dist = shape.Triangle(lower=0.9*baseline, midpoint=baseline, upper=1.1*baseline)
        @param(name='Conditioning temperature', element='HX301', kind='coupled', units='K', baseline=baseline, distribution=dist)
        def set_Conditioning_temperature(i):
            fs_unit.HX301.T = i
        
        dist = shape.Uniform(lower=0.08, upper=1)
        @param(name='Conditioning time', element='CT301', kind='coupled', units='h', baseline=fs_unit.CT301.tau, distribution=dist)
        def set_Conditioning_time(i):
            fs_unit.CT301.tau = i

        dist = shape.Uniform(lower=1, upper=1.5)
        @param(name='Conditioning energy demand', element='CT301', kind='coupled', units='kW/m^3', baseline=fs_unit.CT301.kW_per_m3, distribution=dist)
        def set_Conditioning_energy(i):
            fs_unit.CT301.kW_per_m3 = i

        baseline = fs_unit.DY301.final_moisture
        dist = shape.Triangle(lower=0, midpoint=baseline, upper=0.01)
        @param(name='Fibers moisture', element='DY301', kind='coupled', units='%/100', baseline=baseline, distribution=dist)
        def set_Fibers_moisture(i):
            fs_unit.DY301.final_moisture = i      

        baseline = fs_unit.DY301.drying_T
        dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Drying temeprature', element='DY301', kind='coupled', units='K', baseline=baseline, distribution=dist)
        def set_Drying_temperature(i):
            fs_unit.DY301.drying_T = i    
        
        baseline = fs_unit.DY301.dryer_efficiency
        dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Dryer efficiency', element='DY301', kind='coupled', units='%/100', baseline=baseline, distribution=dist)
        def set_Dryer_efficiency(i):
            fs_unit.DY301.dryer_efficiency = i 

        # --- S4: Module assembly --- (4)

        baseline = fs_unit.MD401.electricity_per_membranearea
        dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Module assembly energy', element='MD401', kind='coupled', units='kWh/m^2 membrane', baseline=baseline, distribution=dist)
        def set_Module_energy(i):
            fs_unit.MD401.electricity_per_membranearea = i 
        
        baseline = fs_unit.MD401.epoxy_per_membranearea
        dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Epoxy resin demand', element='MD401', kind='coupled', units='kg/m^2 membrane', baseline=baseline, distribution=dist)
        def set_Epoxy_resin(i):
            fs_unit.MD401.epoxy_per_membranearea = i 

        if solvent_recycling == 'yes':

            baseline = fs_unit.DT501.tau
            dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
            @param(name='Mixing time, ww', element='DT501', kind='coupled', units='h', baseline=baseline, distribution=dist)
            def set_Mixing_time_ww(i):
                fs_unit.DT501.tau = i 
            
            dist = shape.Uniform(lower=1, upper=1.5)
            @param(name='Mixing energy, ww', element='DT501', kind='coupled', units='kW/m^3', baseline=fs_unit.DT501.kW_per_m3, distribution=dist)
            def set_Mixing_energy_ww(i):
                 fs_unit.DT501.kW_per_m3 = i 

        # --- S5: WW dilution --- (3)

        baseline = fs_unit.WT501.COD_target
        dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        @param(name='Target COD', element='WT501', kind='coupled', units='-', baseline=baseline, distribution=dist)
        def set_Target_COD(i):
            fs_unit.WT501.COD_target = i 
        
        #It does not make sense to test these 2 since this may happen in a concrete infraestructure, so probably no residence time to set (no tanks), or agitation.
        # baseline = fs_unit.WT501.tau                                                      
        # dist = shape.Triangle(lower=baseline*0.9, midpoint=baseline, upper=baseline*1.1)
        # @param(name='Dilution time', element='WT501', kind='coupled', units='h', baseline=baseline, distribution=dist)
        # def set_Dilution_time(i):
        #     fs_unit.WT501.tau = i 
        
        # dist = shape.Uniform(lower=0.01, upper=0.02)
        # @param(name='Dilution energy demand', element='WT501', kind='coupled', units='kW/m^3', baseline=fs_unit.WT501.kW_per_m3, distribution=dist)
        # def set_Dilution_energy_demand(i):
        #     fs_unit.WT501.kW_per_m3 = i 


    def set_lca_params():
        """
        Uncertainty on the background inventories using pedigree approach.
        These are 'isolated', they only change characterization factors.
        It is only included the parameters that show to have an impact bigger than 10%
        in more than 1 indicator, or to to be the major contributor in at leas 1 indicator
        accroding to the baseline analysis.
        """
        # -----------------------------------------------------------------------
        # Using log-normal distribution where:
        #   mu    = 0              (multiplier centered at 1.0 = no change)
        #   sigma = ln(GSD)        
        # Approach: sample a dimensionless multiplier i ~ LogNormal(0, ln(GSD))
        # Then: CF_new = CF_baseline * i  for ALL 18 indicators simultaneously
        # -----------------------------------------------------------------------

        # NMP — all 18 indicators vary together with same GSD²
        NMP_item = qs.ImpactItem.get_item('NMP_item')
        GSD_NMP = 1.056193047   # from pedigree matrix in uncertainty variable.xls
        baseline_CFs_NMP = {ind: NMP_item.CFs[ind] for ind in NMP_item.CFs}  #loops over the indicator key that lives inside qsdsan-the baseline

        dist = shape.LogNormal(mu=0, sigma=np.log(GSD_NMP))   #sigma or variance is the square root of the standar deviation
        @param(name='NMP background', element='LCA', kind='isolated',
            units='-', baseline=1.0, distribution=dist)   #baseline will be the cf_baseline*1 
        def set_NMP_uncertainty(i):
            # i is a multiplier drawn from log-normal centered at 1.0
            # multiply ALL CFs by this factor simultaneously
            for ind, cf_baseline in baseline_CFs_NMP.items():
                NMP_item.CFs[ind] = cf_baseline * i

        # --- Glycerol (conditioning) ---
        glycerol_item         = qs.ImpactItem.get_item('glycerol_item')
        GSD_glycerol          = 1.055238601   
        baseline_CFs_glycerol = {ind: glycerol_item.CFs[ind] for ind in glycerol_item.CFs}

        dist = shape.LogNormal(mu=0, sigma=np.log(GSD_glycerol))
        @param(name='Glycerol background', element='LCA', kind='isolated',
            units='-', baseline=1.0, distribution=dist)
        def set_glycerol_uncertainty(i):
            for ind, cf_baseline in baseline_CFs_glycerol.items():
                glycerol_item.CFs[ind] = cf_baseline * i

        # --- Steam ---
        steam_item         = qs.ImpactItem.get_item('steam_item')
        GSD_steam          = 1.055238601   
        baseline_CFs_steam = {ind: steam_item.CFs[ind] for ind in steam_item.CFs}

        dist = shape.LogNormal(mu=0, sigma=np.log(GSD_steam))
        @param(name='Steam background', element='LCA', kind='isolated',
            units='-', baseline=1.0, distribution=dist)
        def set_steam_uncertainty(i):
            for ind, cf_baseline in baseline_CFs_steam.items():
                steam_item.CFs[ind] = cf_baseline * i

    
        # --- Softened water or water for wastewater dilution---
        water_ww_item         = qs.ImpactItem.get_item('water_wwt_item')
        GSD_water_ww          = 1.060597437   
        baseline_CFs_waterww = {ind: water_ww_item.CFs[ind] for ind in water_ww_item.CFs}

        dist = shape.LogNormal(mu=0, sigma=np.log(GSD_water_ww))
        @param(name='Water in wastewater background', element='LCA', kind='isolated',
            units='-', baseline=1.0, distribution=dist)
        def set_water_ww_uncertainty(i):
            for ind, cf_baseline in baseline_CFs_waterww.items():
                water_ww_item.CFs[ind] = cf_baseline * i
       
        # --- Wastewater treatment burden to the wastewater produced ---
        wastewater_item         = qs.ImpactItem.get_item('wastewater_item')
        GSD_wastewater          = 1.06596836  # replace with your pedigree value
        baseline_CFs_wastewater = {ind: wastewater_item.CFs[ind] for ind in wastewater_item.CFs}

        dist = shape.LogNormal(mu=0, sigma=np.log(GSD_wastewater))
        @param(name='Wastewater treatment background', element='LCA', kind='isolated',
            units='-', baseline=1.0, distribution=dist)
        def set_wastewater_uncertainty(i):
            for ind, cf_baseline in baseline_CFs_wastewater.items():
                wastewater_item.CFs[ind] = cf_baseline * i
       
        if polymer_option in ('CA', 'CA_bioAA'):
            # --- Ethanol ---
            ethanol_item         = qs.ImpactItem.get_item('ethanol_item')
            GSD_ethanol          = 1.056193047   
            baseline_CFs_ethanol = {ind: ethanol_item.CFs[ind] for ind in ethanol_item.CFs}

            dist = shape.LogNormal(mu=0, sigma=np.log(GSD_ethanol))
            @param(name='Ethanol background', element='LCA', kind='isolated',
                units='-', baseline=1.0, distribution=dist)
            def set_ethanol_uncertainty(i):
                for ind, cf_baseline in baseline_CFs_ethanol.items():
                    ethanol_item.CFs[ind] = cf_baseline * i
    

        # Load fitted parameters
        fitted = pd.read_csv(
            os.path.join(os.path.dirname(__file__), 'fitted_params_all_polymers.csv')
        )

        # =========================================================================
        # PART 2 — Polymer CF uncertainty
        # =========================================================================

        # --- PSF module housing — ALWAYS included regardless of polymer_option ---
        # Because module housing is always PSF even when membrane is CA or CA_bioAA
        PSF_module_item         = qs.ImpactItem.get_item('PSF_module_item')

        for ind in INDICATOR_IDS:
            row = fitted[(fitted['polymer'] == 'PSF') &
                        (fitted['indicator'] == ind)]
            if len(row) == 0:
                continue

            mu_fit    = float(row['mu'].values[0])
            sigma_fit = float(row['sigma'].values[0])

            dist = shape.LogNormal(mu=mu_fit, sigma=sigma_fit)

            @param(name=f'PSF module CF {ind}', element='LCA', kind='isolated',
                units='kg CF/kg polymer',
                baseline=float(np.exp(mu_fit)),
                distribution=dist)
            def make_setter_PSF_module(ind=ind):
                def set_PSF_module_CF(i):
                    PSF_module_item.CFs[ind] = i
                return set_PSF_module_CF
            make_setter_PSF_module(ind)

        # --- Membrane polymer — depends on polymer_option ---

        if polymer_option == 'PSF':
            PSF_item         = qs.ImpactItem.get_item('PSF_item')

            for ind in INDICATOR_IDS:
                row = fitted[(fitted['polymer'] == 'PSF') &
                            (fitted['indicator'] == ind)]
                if len(row) == 0:
                    continue

                mu_fit    = float(row['mu'].values[0])
                sigma_fit = float(row['sigma'].values[0])

                dist = shape.LogNormal(mu=mu_fit, sigma=sigma_fit)

                @param(name=f'PSF membrane CF {ind}', element='LCA', kind='isolated',
                    units='kg CF/kg polymer',
                    baseline=float(np.exp(mu_fit)),
                    distribution=dist)
                def make_setter_PSF(ind=ind):
                    def set_PSF_CF(i):
                        PSF_item.CFs[ind] = i
                    return set_PSF_CF
                make_setter_PSF(ind)

        elif polymer_option == 'CA':
            CA_item         = qs.ImpactItem.get_item('CA_item')

            for ind in INDICATOR_IDS:
                row = fitted[(fitted['polymer'] == 'CA') &
                            (fitted['indicator'] == ind)]
                if len(row) == 0:
                    continue

                c_fit     = float(row['shape'].values[0])
                loc_fit   = float(row['loc'].values[0])
                scale_fit = float(row['scale'].values[0])

                dist = shape.Weibull(c_fit, scale=scale_fit) + loc_fit

                @param(name=f'CA membrane CF {ind}', element='LCA', kind='isolated',
                    units='kg CF/kg polymer',
                    baseline=float(row['median'].values[0]),
                    distribution=dist)
                def make_setter_CA(ind=ind):
                    def set_CA_CF(i):
                        CA_item.CFs[ind] = i
                    return set_CA_CF
                make_setter_CA(ind)

        elif polymer_option == 'CA_bioAA':
            CA_bioAA_item         = qs.ImpactItem.get_item('CA_bioAA_item')

            for ind in INDICATOR_IDS:
                row = fitted[(fitted['polymer'] == 'CA_bioAA') &
                            (fitted['indicator'] == ind)]
                if len(row) == 0:
                    continue

                a_fit     = float(row['shape'].values[0])
                loc_fit   = float(row['loc'].values[0])
                scale_fit = float(row['scale'].values[0])

                dist = shape.Gamma(a_fit, scale=scale_fit) + loc_fit

                @param(name=f'CA_bioAA membrane CF {ind}', element='LCA',
                    kind='isolated',
                    units='kg CF/kg polymer',
                    baseline=float(row['median'].values[0]),
                    distribution=dist)
                def make_setter_bioAA(ind=ind):
                    def set_bioAA_CF(i):
                        CA_bioAA_item.CFs[ind] = i
                    return set_bioAA_CF
                make_setter_bioAA(ind)

    if analysis == 'all':
        set_technological_params()
        set_lca_params()
    elif analysis == 'technological':
        set_technological_params()
    elif analysis == 'background':
        set_lca_params()
    else:
        raise RuntimeError(f'In create_model(sys,analysis,membrane_area_per_year, polymer_option, solvent_recycling),parameter={analysis} is not "all", "technological" or "background". Please define as one of these.')
        

    # =========================================================================
    # METRICS
    # =========================================================================
    # Metrics tell the model WHAT TO MEASURE after each simulation run.
    # Each metric needs:
    #   name    : readable label for plots and tables
    #   units   : string for axis labels
    #   element : 'LCA' 
    #
    # -------------------------------------------------------------------------

    @metric(name='Global Warming', units='kg CO2-Eq/m²', element='LCA')
    def get_GWP():
        return lca.get_total_impacts(annual=True)['GWP'] / membrane_area_per_year

    @metric(name='Terrestrial Acidification', units='kg SO2-Eq/m²', element='LCA')
    def get_TAP():
        return lca.get_total_impacts(annual=True)['TAP'] / membrane_area_per_year

    @metric(name='Freshwater Ecotoxicity', units='kg 1,4-DCB-Eq/m²', element='LCA')
    def get_FETP():
        return lca.get_total_impacts(annual=True)['FETP'] / membrane_area_per_year

    @metric(name='Marine Ecotoxicity', units='kg 1,4-DCB-Eq/m²', element='LCA')
    def get_METP():
        return lca.get_total_impacts(annual=True)['METP'] / membrane_area_per_year

    @metric(name='Terrestrial Ecotoxicity', units='kg 1,4-DCB-Eq/m²', element='LCA')
    def get_TETP():
        return lca.get_total_impacts(annual=True)['TETP'] / membrane_area_per_year

    @metric(name='Fossil Fuel Potential', units='kg oil-Eq/m²', element='LCA')
    def get_FFP():
        return lca.get_total_impacts(annual=True)['FFP'] / membrane_area_per_year

    @metric(name='Freshwater Eutrophication', units='kg P-Eq/m²', element='LCA')
    def get_FEP():
        return lca.get_total_impacts(annual=True)['FEP'] / membrane_area_per_year

    @metric(name='Marine Eutrophication', units='kg N-Eq/m²', element='LCA')
    def get_MEP():
        return lca.get_total_impacts(annual=True)['MEP'] / membrane_area_per_year

    @metric(name='Human Toxicity Carcinogenic', units='kg 1,4-DCB-Eq/m²', element='LCA')
    def get_HTPc():
        return lca.get_total_impacts(annual=True)['HTPc'] / membrane_area_per_year

    @metric(name='Human Toxicity Non-carcinogenic', units='kg 1,4-DCB-Eq/m²', element='LCA')
    def get_HTPnc():
        return lca.get_total_impacts(annual=True)['HTPnc'] / membrane_area_per_year

    @metric(name='Ionising Radiation', units='kBq Co-60-Eq/m²', element='LCA')
    def get_IRP():
        return lca.get_total_impacts(annual=True)['IRP'] / membrane_area_per_year

    @metric(name='Land Use', units='m²a crop-Eq/m²', element='LCA')
    def get_LOP():
        return lca.get_total_impacts(annual=True)['LOP'] / membrane_area_per_year

    @metric(name='Surplus Ore Potential', units='kg Cu-Eq/m²', element='LCA')
    def get_SOP():
        return lca.get_total_impacts(annual=True)['SOP'] / membrane_area_per_year

    @metric(name='Ozone Depletion', units='kg CFC-11-Eq/m²', element='LCA')
    def get_ODP():
        return lca.get_total_impacts(annual=True)['ODP'] / membrane_area_per_year

    @metric(name='Particulate Matter Formation', units='kg PM2.5-Eq/m²', element='LCA')
    def get_PMFP():
        return lca.get_total_impacts(annual=True)['PMFP'] / membrane_area_per_year

    @metric(name='Photochemical Ox. Human Health', units='kg NOx-Eq/m²', element='LCA')
    def get_HOFP():
        return lca.get_total_impacts(annual=True)['HOFP'] / membrane_area_per_year

    @metric(name='Photochemical Ox. Ecosystems', units='kg NOx-Eq/m²', element='LCA')
    def get_EOFP():
        return lca.get_total_impacts(annual=True)['EOFP'] / membrane_area_per_year

    @metric(name='Water Consumption', units='m³/m²', element='LCA')
    def get_WCP():
        return lca.get_total_impacts(annual=True)['WCP'] / membrane_area_per_year

    return model