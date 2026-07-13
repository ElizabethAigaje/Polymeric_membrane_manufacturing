"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 04/27/2026
"""
# %%
#Import packages

import qsdsan as qs
import biosteam as bst
from biosteam.units import HXutility
import numpy as np
import pandas as pd
import os
from logging import raiseExceptions
import matplotlib.pyplot as plt

#import functions from other files
from components_membrane import create_components_membrane
from _lca_data_membrane_ import _load_lca_data_membrane
from _lca_results_ import LCAResults
from _model_membrane_ import *
from _uncertainty_analysis_ import uncertainty_analysis
from _sensitivity_analysis_ import sensitivity_analysis

#import unit operations
from Units.mixer_solvent import Mixer_solvent
from Units.mixer_additive import Mixer_additive
from Units.mixer_water import Mixer_water
from Units.Polymersol_Tank import Dissolution_tank
from Units.Polymersol_DegassTank import Dissolution_degastank
from Units.Boresol_DegassTank import Bore_degasTank
from Units.Extrusion import Extruder
from Units.mixer_non_solvent import Mixer_non_solvent
from Units.Coagulation import Coagulation_bath
from Units.mixer_rinsing_water import Mixer_rinsing_water
from Units.Rinsing_tank import Rinsing_tank
from Units.mixer_glycerol_conditioning import Mixer_glycerol_conditioning
from Units.mixer_water_conditioning import Mixer_water_conditioning
from Units.Preconditioning_tank import Conditioning_tank
from Units.Dryer import Dryer
from Units.Module_assembly import Module_assembly
from Units.mixer_solventr import Mixer_solvent_r
from Units.Regeneration_Tank import Regeneration_tank
from Units.Regeneration_bath import Regeneration_bath
from Units.Dilution_Tank import Dilution_tank

#For LCA
cf_dct = _load_lca_data_membrane()


os.environ["PATH"] += os.pathsep + r'C:\Users\eka5489\AppData\Local\anaconda3\envs\REEPS_env\Library\bin'
bst.nbtutorial()

#Collect components
cmps1= create_components_membrane()

def system_membrane(polymer_solution, polymer_option, solvent_recycling):  
    """ Creates the class that models the membrane manufacturing process and 
        performs a life cycle assessment (LCA) using Ecoinvent 3.10 database for background 
        data and the ReCiPe 2016 v1.03 midpoint H cut-off method.

    Parameters
    ----------
    polymersol : float
        Volume of polymer solution processed per day, [L/day] (the default is 2000 L).
    polumer_option : str
        PSF, CA, CA_bioAA, referring to the polysulfone, cellulose acetate, cellulose acetate prepared using biobased acetic acid, respectively.
    solvent_recycling : str
        Wheter or not the system includes solvent recovery and recycling, the options are: yes or no.
    
    Returns
    -----
    system : qsdsan.System
        Process simulated, mass flows in and out of the system, process diagram.
    lca : qsdsan.lca
        Table of given annual impact categories per construction, transport, streams and others.
    membrane_area_per_year : float
        Annual membrane output [m²/year], used as the functional unit
        for normalizing LCA results. Calculated from polymer
        mass flow rate and polymer mass per membrane area [kg/m²].

    See Also
    --------
    class MembraneParams: Class with its description.
        
    """

    #-----SYTEM COMPONENTS ----

    polymer_database = {
    'PSF': {
        'component': cmps1['polysulfone'],
        'density': 1240,  # kg/m3
        'post_treatment': None,
    },

    'CA': {
        'component': cmps1['polysulfone'],  #we are gonna treat all the other properties same as PSF, except for the density
        'density': 1300,  # kg/m^3 (range 1280 to 1320 kg/m^3)
        'post_treatment': 'regeneration',  
    },

    'CA_bioAA': {
        'component': cmps1['polysulfone'],  #we are gonna treat all the other properties same as PSF, except for the density
        'density': 1300,  # kg/m^3 (range 1280 to 1320 kg/m^3)
        'post_treatment': 'regeneration',  
    }
}
    
    if polymer_option not in polymer_database:
        raise ValueError(f"polymer_option must be one of {list(polymer_database.keys())}")

    polymer_data = polymer_database[polymer_option]
    d_polymer = polymer_data['density']
    

    solvent=cmps1['NMP']
    additive=cmps1['PEG']
    water = cmps1['H2O']
    glycerol =cmps1['glycerol']

    #-----MEMBRANE SYSTEM PARAMETERS ----

    class MembraneParams:
        """
        Parameters for the membrane manufacturing system.

        This class stores all system-level variables used across unit operations,
        including independent inputs, operating conditions,
        and derived properties that are computed based on these inputs.

        Parameters are updated dynamically through the `update()` method to ensure
        consistency across all units during simulation and uncertainty analysis.
        ...

        Attributes
        ----------
        
        --- Independent variables ---
        polymer_fraction : float
            Mass fraction of polymer in the polymer solution, [w/w %/100] (the default is 0.18 [1]).
        additive_fraction : float
            Mass fraction of solvent in the polymer solution, [w/w %/100] (the default is 0.10 [1-4]).
        v_polymersol : float
            Volume of polymer solution processed per day, [L/day] (Tanks of 2000 L considered).

        --- Operating conditions ---
        T_polymersol : float
            Temperature of polymer solution preparation, [°C] (the default is 60°C [5,6]).
        T_boresol : float
            Temperature of bore solution preparation, [°C] (the default is 60°C).
        
        --- Material and system constants ---
        d_polymer : float
            Density of polymer [kg/m³].
        ratio_boresol_to_polymersol : float
            Volume ratio of bore solution to polymer solution processed per day, [-] (the default is 1).
        solventbore_fraction : float
            Mass fraction of solvent in the bore solution, [-] (the default is 0.5  [6]).
        conditioning_concentration : float
            Glycerol mass concentration in the aqueous solution, [w/w% /100], (the default is 0.25*).

        --- Fiber/module geometry (based on references 7, 8, 9)---
        Di : float
            Inner fiber diameter, [mm] (the default is 0.93).
        Do : float
            Outer fiber diameter, [mm] (the default is 1.67).
         L : float
            Fiber length, [m] (the default is 1.2).
        porosity : float
            Membrane porosity, [-] (the default is 0.7*).
        fibernumber_per_module : int
            Number of fibers per module, [-] (the default is 2600).
        membranearea_per_module : float
            Membrane area per module,[m²] (the default is 9).
        Do_module : float
            Module outer diameter [m] (the default is 0.18).
        L_module : float
            Module length [m] (the default is 1.2).
        thickness_module : float
            Module wall thickness [m] (the default is 0.005).

        --- Derived properties (computed in update) ---
        d_polymersol : float
            Density of polymer solution [kg/m³].
        d_boresol : float
            Density of bore solution [kg/m³].
        ratio_solvent : float
            Fraction of solvent allocated to polymer solution [-].
        d_conditioningsol : float
            Density of conditioning solution [kg/m³].
        polymermass_to_membranearea : float
            Polymer mass per membrane area [kg/m²].
        cartridge_mass_per_area : float
            Cartridge polymer material per membrane area [kg/m²].

        Methods
        -------
        update(solvent, additive, water, glycerol)
            Recalculate all dependent properties based on current parameter values.
        
        Notes: 
        ------
        *Based on the average range provided by the industry members of the team.

        References
        ----------
        1. Q. Zheng, et al. The relationship between porosity and kinetics parameter of  membrane formation in PSF ultrafiltration membrane. Journal of Membrane Science 286 (2006) 7–11
        2. W. Zhao, et al. Facile method to fabricate robust homogeneous braid‐reinforced cellulose acetate hollow fiber membranes with enhanced physiochemical properties. 
           Cellulose (2024) 31:395–410 
        3. X. Gao, et al. Fabrication and characterization of cellulose acetate  ultrafiltration membrane and its application for efficient  bovine serum albumin separation. 
           Polym Eng Sci. 2023;63:1423–1438.
        4. S. Zheng, et al.Performance investigation of hydrophilic regenerated cellulose ultrafiltration membranes with excellent anti-fouling property via hydrolysis technology. 
           Journal of Environmental Chemical Engineering 11 (2023) 109041.
        5. R. Chidambaran, et al. United States Patent Application Publication, US 2012/0111790a1, 2012.
        6. Y. Bang, et al. Influence of bore fluid composition on the physiochemical properties and performance of hollow fiber membranes for ultrafiltration. Chemosphere 259 (2020) 127467
        7. Cytiva. MaxCellTM ultrafiltration hollow fiber tangential flow cartridges: UFP-300-E-85. Cytiva.com
        8. Cytiva. MaxCell process-scale hollow fiber cartridges: Specifications. Cytiva.com
        9. F. Prezelus, et al. A generic process modelling – LCA approach for UF membrane fabrication: Application to cellulose acetate membranes. Journal of Membrane Science 618 (2021) 118594

        """
        def __init__(self):

            # -- independent (uncertainty)
            self.polymer_fraction = 0.18   
            self.additive_fraction = 0.10   
            self.v_polymersol = polymer_solution

            #-- constants --
            self.d_polymer = d_polymer
            self.T_polymersol = 60  
            self.T_boresol = 60     

            self.ratio_boresol_to_polymersol = 1   
            self.solventbore_fraction = 0.5      

            self.conditioning_concentration = 0.25
            
            self.Di = 0.93
            self.Do = 1.67
            self.L = 1.2
            self.porosity = 0.7
            self.fibernumber_per_module = 2600
            self.membranearea_per_module = 9
            self.Do_module = 0.18  #m
            self.L_module= 1.2 #m
            self.thickness_module = 0.005 #m
            self.volume_per_fiber = (3.1416/4)*((self.Do**2)-(self.Di**2))*self.L/(1000**2)

            #-- placeholders (will be computed)
            self.d_polymersol = None
            self.d_boresol = None
            self.ratio_solvent = None
            self.d_conditioningsol = None
            self.polymermass_to_membranearea = None
            self.cartridge_mass_per_area = None

            # --- placeholders for refresh() 
            self._solvent  = None
            self._additive = None
            self._water    = None
            self._glycerol = None


        def update(self, solvent, additive, water, glycerol):
            
            # --- store references so refresh() can call update() later 
            self._solvent = solvent
            self._additive = additive
            self._water    = water
            self._glycerol = glycerol
            
            polymer_fraction = self.polymer_fraction
            additive_fraction = self.additive_fraction
            solvent_fraction = 1 - polymer_fraction - additive_fraction

        # --- polymer solution density ---
            self.d_polymersol = (
                solvent.rho(phase='l', T=self.T_polymersol+273.15, P=101325) * solvent_fraction +
                additive.rho(phase='l', T=self.T_polymersol+273.15, P=101325) * additive_fraction +
                self.d_polymer * polymer_fraction
            )

        # --- bore solution density ---
            water_fraction = 1 - self.solventbore_fraction
            self.d_boresol = (
                solvent.rho(phase='l', T=self.T_boresol+273.15, P=101325) * self.solventbore_fraction +
                water.rho(phase='l', T=self.T_boresol+273.15, P=101325) * water_fraction
            )

        # --- solvent ratio ---
            self.ratio_solvent = (
                self.v_polymersol * self.d_polymersol * solvent_fraction) / (
                self.v_polymersol * self.ratio_boresol_to_polymersol * self.d_boresol * self.solventbore_fraction
                + self.v_polymersol * self.d_polymersol * solvent_fraction
            )

        # --- fiber geometry ---

            self.polymermass_per_fiber = self.volume_per_fiber * self.d_polymer * (1 - self.porosity)
            self.polymermass_to_membranearea = (
                self.polymermass_per_fiber * self.fibernumber_per_module / self.membranearea_per_module
            )

        #----conditioning

            self.d_conditioningsol=(
                glycerol.rho(phase='l', T=25+273.5, P=101324)*(self.conditioning_concentration) +
                water.rho(phase='l', T=25+273.5, P=101324)*(1-self.conditioning_concentration))

        # --- cartridge  ---
            self.cartridge_mass_per_area = ((3.1416*self.L_module*self.d_polymer*((self.Do_module**2)-(self.Do_module-self.thickness_module)**2))/4/
                                            self.membranearea_per_module)  #kg/m^2, densitiy of PSF =1240
        def refresh(self):

            """Re-run update() using stored component references.
            Call this after changing any parameter attribute from outside."""

            self.update(self._solvent, self._additive, 
                        self._water,  self._glycerol)
   
    params = MembraneParams()
    params.update(solvent, additive, water, glycerol)  # call once at the system creation to calcualte derived properties and store _solvent, ....

    
    #-----CREATING INLET STREAMS-----
    # Conditional in membrane polymer option
    if polymer_option == 'PSF':
        qs.SanStream('polymer', stream_impact_item=qs.StreamImpactItem.get_item('PSF_item'), price=15, phase='s', T=25+273, P=101325)  
    elif polymer_option == 'CA':
        qs.SanStream('polymer', stream_impact_item=qs.StreamImpactItem.get_item('CA_item'), price=3.873, phase='s', T=25+273, P=101325)  #Source https://www.intratec.us/solutions/primary-commodity-prices/commodity/cellulose-acetate-prices
    elif polymer_option == 'CA_bioAA':
        qs.SanStream('polymer', stream_impact_item=qs.StreamImpactItem.get_item('CA_bioAA_item'), price=3.873, phase='s', T=25+273, P=101325)

    qs.SanStream('solvent', stream_impact_item=qs.StreamImpactItem.get_item('NMP_item'), price=4.07, phase='l', T=25+273, P=101325)  #$/kg
    qs.SanStream('additive', stream_impact_item=qs.StreamImpactItem.get_item('PEG_item'), price=1.43, phase='l', T=25+273, P=101325)  
    qs.SanStream('nitrogen',stream_impact_item=qs.StreamImpactItem.get_item('nitrogen_item'), price=0.3758, phase='g', T=25+273, P=100000) #Source https://puritygas.ca/revisiting-the-costs-of-nitrogen-gas/#:~:text=According%20to%20Purity%20Gas%2C%20the%20cost%20of,*%20Federally%20imposed%20fossil%20fuel%20carbon%20taxes
    qs.SanStream('glycerol_conditioning', stream_impact_item=qs.StreamImpactItem.get_item('glycerol_item'), price=4.125, phase='l', T=25+273, P=101325) 
    qs.SanStream('polysulfone_module', stream_impact_item=qs.StreamImpactItem.get_item('PSF_module_item'), price=15, phase='l', T=25+273, P=101325) 
    qs.SanStream('epoxy_module', stream_impact_item=qs.StreamImpactItem.get_item('epoxy_item'), price=3.685, phase='l', T=25+273, P=101325) 
    
    # DI and tap water streams
    qs.SanStream('water_boresol', stream_impact_item=qs.StreamImpactItem.get_item('water_boresol_item'), price= 0.001103,  phase='l', T=25+273, P=101325) 
    qs.SanStream('non_solvent', stream_impact_item=qs.StreamImpactItem.get_item('nonsolvent_item'), price= 0.001103, phase='l',T=25+273, P=101325)
    qs.SanStream('rinsing_water',stream_impact_item=qs.StreamImpactItem.get_item('water_rinsing_item'), price= 0.001103, phase='l',T=25+273, P=101325)
    qs.SanStream('water_conditioning', stream_impact_item=qs.StreamImpactItem.get_item('water_conditioning_item'), price= 0.001103, phase='l', T=25+273, P=101325) 
    
    # CA-only streams
    qs.SanStream('sodium_hydroxide', stream_impact_item=qs.StreamImpactItem.get_item('NaOH_item'), price= 0.36, phase='s', T=25+273, P=101325)   #Source https://businessanalytiq.com/procurementanalytics/index/sodium-hydroxide-price-index/
    qs.SanStream('ethanol_rsolution', stream_impact_item=qs.StreamImpactItem.get_item('ethanol_item'), price= 0.642, phase='l', T=25+273, P=101325)    #Source https://tradingeconomics.com/commodity/ethanol#:~:text=Ethanol%20rose%20to%201.92%20USD/Gal%20on%20April,news%20%2D%20updated%20on%20April%20of%202026.
    qs.SanStream('water_rsolution', stream_impact_item=qs.StreamImpactItem.get_item('water_rsolution_item'), price= 0.001103, phase='l', T=25+273, P=101325)
    qs.SanStream('water_wwt', stream_impact_item=qs.StreamImpactItem.get_item('water_wwt_item'), price= 0.00037, phase='l', T=25+273, P=101325) 

    #-----CREATING WASTE STREAMS-----
    #Note: negative CF in the item already hanldes the sign correctly (as waste)
    qs.WasteStream('diluted_wastewater', stream_impact_item=qs.StreamImpactItem.get_item('wastewater_item'), price= 0.000056, phase='l', T=25+273, P=101325) # $/kg in Turton, R, et al., Analysis, synthesis and design of chemical processes, pag.246 (56$/1000m^3)

    # --------------------------------------------------------------------------
    # CREATE SYSTEM
    # --------------------------------------------------------------------------

    #use qs.Flowsheet serves as a centrilized structure that manages all unit operations,streams and system object in a process system.
    flowsheet = qs.Flowsheet.flowsheet.default 
    fs_stream = flowsheet.stream 
    fs_unit = flowsheet.unit 

    #Polymer solution
    if solvent_recycling == 'no':

        #Polymer solution
        M1=Mixer_solvent('M101', ins=fs_stream.solvent, outs='solvent101', params =params, rigorous=False)
        ST1=qs.sanunits.StorageTank('ST101', ins= M1-0, outs='solvent_c', tau=336, 
                                vessel_material='Stainless steel', vessel_type='Cone roof')  #low volatility
        P1=qs.sanunits._pumping.Pump('P101', ins=ST1-0, outs='solvent102', P=101324, material='Stainless steel')
        HX1=qs.sanunits.HXutility('HX101', ins=P1-0, outs='solvent_h', T= params.T_polymersol+273.15, rigorous=False)
        S1=qs.sanunits.Splitter('S101', ins=HX1-0, outs=('solvent_polymersol', 'solvent_bosol'), split=params.ratio_solvent)  
       
       
        M2=Mixer_additive('M102', ins=fs_stream.additive, outs='additive101', params=params, rigorous=False)
        ST2=qs.sanunits.StorageTank('ST102', ins= M2-0, outs='additive_c', tau=336, 
                                    vessel_material='Stainless steel', vessel_type='Cone roof')  
        P2=qs.sanunits._pumping.Pump('P102', ins=ST2-0, outs='additive102', P=101324, material='Stainless steel')
        HX2=qs.sanunits.HXutility('HX102', ins=P2-0, outs='additive_h', T= params.T_polymersol+273.15, rigorous=False)  
    
        P3=qs.sanunits._pumping.Pump('P103', ins=S1-0, outs='solvent103', P=101324, material='Stainless steel')
        P4=qs.sanunits._pumping.Pump('P104', ins=HX2-0, outs='additive103', P=101324, material='Stainless steel')
    
        DT1= Dissolution_tank('DT101', ins=(fs_stream.polymer, P3-0, P4-0), outs=('polymer_solution'), params=params, tau=24, kW_per_m3=14)
        P5=qs.sanunits._pumping.Pump('P105', ins=DT1-0, outs='polymer_ndegassolution', P=101324, material='Stainless steel')
        DG1 = Dissolution_degastank('DG101', ins=P5-0, outs=('polymer_degassolution'), tau=12, kW_per_m3=0.07, P=700)   

        #Bore solution
        M3= Mixer_water('M103', ins=fs_stream.water_boresol, outs='water_c', params=params, rigorous=False)
        HX3=qs.sanunits.HXutility('HX103', ins=M3-0, outs='water_h', T= params.T_boresol+273.15, rigorous=False)
        P6=qs.sanunits._pumping.Pump('P106', ins=HX3-0, outs='water101', P=101324, material='Stainless steel')
    
        P7=qs.sanunits._pumping.Pump('P107', ins=S1-1, outs='solvent104', P=101324, material='Stainless steel')
        DT2 = qs.sanunits.MixTank('DT102', ins=( P6-0, P7-0), outs=('bore_solution'), tau=1, kW_per_m3=1.25)  #Ranges from 1.0-1.5, which correspond to liquid-liquid mixing 
        P8=qs.sanunits._pumping.Pump('P108', ins=DT2-0, outs='bore_ndegassolution', P=101324, material='Stainless steel')
        DG2 = Bore_degasTank ('DG102', ins=(P8-0, fs_stream.nitrogen), outs=('bore_degassolution', 'nitrogen_out'), nitrogen_demand= 0.00056, params=params, tau=1, kW_per_m3=0.07)
        
        #Extrusion, coagulation & rinsing
        EX1= Extruder('EX201', ins=(DG1-0, DG2-0), outs=('lumen'), params=params, power_demand = 300)
        M4= Mixer_non_solvent('M201', ins=fs_stream.non_solvent, outs='non_solvent201', nonsolvent_per_membranearea=150, params=params, T=25+273, rigorous=False)
        HX4=qs.sanunits.HXutility('HX201', ins=M4-0, outs='nonsolvent_h', T= 20+273.15, rigorous=False)
        CB1= Coagulation_bath('CB201', ins=(EX1-0, HX4-0), outs=('precipitated_lumen', 'wastewater1'), solvent_out=0.7)
        
        if polymer_option in ('CA', 'CA_bioAA'):
            RGT1= Regeneration_tank('RGT201', ins=(fs_stream.sodium_hydroxide, fs_stream.ethanol_rsolution, fs_stream.water_rsolution,'ethanolsol_recycled'), outs='regeneration_solution',
                 params=params, naoh_concentration=0.2, tau=0.5, kW_per_m3=2)   
            #HXB4= qs.sanunits.HXutility('HX201b', ins=RGT1-0, outs='regeneration_solutionh', T= 70+273.15, rigorous=False) #seems that regeneration is done at ambient temperature
            RB1= Regeneration_bath('RB201', ins=(CB1-0, RGT1-0), outs=('lumen_for_rinsing', 'wastewater1b'), tau=0.5, kW_per_m3=1.25) #from Ref 1 in regeneration tank
            SP1 = qs.sanunits.Splitter ('SP201', ins=RB1-1, outs=('purge', 3-RGT1), split=0.025)  #use a purge so Na-aacetate does not build up in the stream that is concentrated with EtOH, well controled systems purge 1-5%

            stream_to_rinse = RB1-0
        else:
            stream_to_rinse = CB1-0
        
        M5= Mixer_rinsing_water('M202', ins=fs_stream.rinsing_water, outs='rinsing_water201', rinsing_per_membranearea=150, params=params, rigorous=False)
        HX5=qs.sanunits.HXutility('HX202', ins=M5-0, outs='rinsing_waterh', T= 50+273.15, rigorous=False)
        RT1= Rinsing_tank('RT201', ins=(stream_to_rinse, HX5-0), outs= ('lumen_after_rinsing', 'wastewater2'), solvent_out=0.99999, tau=0.5)   #test 30 min range 5 min to 1 h
        
        #Preconditioning and drying
        M6= Mixer_glycerol_conditioning('M301', ins=fs_stream.glycerol_conditioning, outs='glycerol_conditioning301',  
                                        fraction_filled=1.1, params=params, rigorous=False)
        ST3=qs.sanunits.StorageTank('ST301', ins= M6-0, outs='glycerol_conditioning302', tau=336, 
                                    vessel_material='Stainless steel', vessel_type='Cone roof')  #low volatility
        P9=qs.sanunits._pumping.Pump('P301', ins=ST3-0, outs='glycerol_conditioning303', P=101324, material='Stainless steel')  #to the tank

        M7= Mixer_water_conditioning('M302', ins=fs_stream.water_conditioning, outs='water_conditioning301',  
                                    fraction_filled=1.1, params=params, rigorous=False)
        P10=qs.sanunits._pumping.Pump('P302', ins=M7-0, outs='water_conditioning301', P=101324, material='Stainless steel')  #to the tank

        DT3 = qs.sanunits.MixTank('MT301', ins=( P9-0, P10-0), outs=('conditioning_solutionc'), tau=0.5, kW_per_m3=1.25)  #Ranges from 1.0-1.5, which correspond to liquid-liquid mixing ##Residence time suggested for chta gpt for easy mixtures
        HX6=qs.sanunits.HXutility('HX301', ins=DT3-0, outs='conditioning_solutionh', T = 30+273, rigorous=False)
        P11=qs.sanunits._pumping.Pump('P303', ins=HX6-0, outs='conditioning_solution301', P=101324, material='Stainless steel')
        CT1 = Conditioning_tank('CT301', ins=(RT1-0, P11-0), outs=('conditioned_fiber'),tau=0.5, kW_per_m3=1.25, conditioning_temperature=30+273)   #same as rinsing tank
        
        DY1 = Dryer('DY301', ins=(CT1-0), outs=('solid_fiber', 'wastewater3'), final_moisture=0.0, drying_T = 50+273.15, dryer_efficiency = 0.7)
        HX7 = qs.sanunits.HXutility('HX302', ins=DY1-1, outs='wastewater3', T= 25+273.15, rigorous=False)
        
        #Module assembly
        MD1 = Module_assembly('MD401', ins=(DY1-0, fs_stream.polysulfone_module, fs_stream.epoxy_module), outs=('module'), electricity_per_membranearea= 0.0667, 
                            epoxy_per_membranearea=0.127, params=params)
        
        #Wastewater treatment
        if polymer_option in ('CA', 'CA_bioAA'): 
            M8 = qs.sanunits.Mixer('M501', ins =(CB1-1, RT1-1, SP1-0, HX7-0), outs='wastewater')    
        else:
            M8 = qs.sanunits.Mixer('M501', ins =(CB1-1, RT1-1, HX7-0), outs='wastewater')
        
        WT1 = Dilution_tank('WT501', ins=(M8-0, fs_stream.water_wwt), outs=fs_stream.diluted_wastewater, COD_target=425.34, tau=0.5, kW_per_m3=0.0)  

        if polymer_option in ('CA', 'CA_bioAA'):
            sys_membrane=qs.System('sys_membrane', path=(M1, ST1, P1, HX1, S1, M2, ST2, P2, HX2, P3, P4, DT1, P5, DG1, M3, HX3, P6, P7, DT2, P8, DG2, EX1,M4, HX4, CB1, RGT1, 
                                                         RB1, SP1, M5, HX5, RT1, M6, M7, P9, ST3, P10, DT3, HX6, P11, CT1, DY1, HX7, MD1, M8, WT1), recycle = SP1-1)
        else:
            sys_membrane=qs.System('sys_membrane', path=(M1, ST1, P1, HX1, S1, M2, ST2, P2, HX2, P3, P4, DT1, P5, DG1, M3, HX3, P6, P7, DT2, P8, DG2, EX1,M4, HX4, CB1, M5, 
                                                         HX5, RT1, M6, M7, P9, ST3, P10, DT3, HX6, P11, CT1, DY1, HX7, MD1, M8,WT1))
       
       
    elif solvent_recycling == 'yes':

        #Polymer solution
        M1=Mixer_solvent_r('M101', ins=(fs_stream.solvent, 'solvent_recycle'), outs='solvent101', params=params, rigorous=False)
        ST1=qs.sanunits.StorageTank('ST101', ins= M1-0, outs='solvent_c', tau=336, 
                                vessel_material='Stainless steel', vessel_type='Cone roof')  #low volatility
        P1=qs.sanunits._pumping.Pump('P101', ins=ST1-0, outs='solvent102', P=101324, material='Stainless steel')
        HX1=qs.sanunits.HXutility('HX101', ins=P1-0, outs='solvent_h', T= params.T_polymersol+273.15, rigorous=False)
        S1=qs.sanunits.Splitter('S101', ins=HX1-0, outs=('solvent_polymersol', 'solvent_bosol'), split=params.ratio_solvent)

        M2=Mixer_additive('M102', ins=fs_stream.additive, outs='additive101', params=params, rigorous=False)
        ST2=qs.sanunits.StorageTank('ST102', ins= M2-0, outs='additive_c', tau=336, 
                                    vessel_material='Stainless steel', vessel_type='Cone roof')  
        P2=qs.sanunits._pumping.Pump('P102', ins=ST2-0, outs='additive102', P=101324, material='Stainless steel')
        HX2=qs.sanunits.HXutility('HX102', ins=P2-0, outs='additive_h', T= params.T_polymersol+273.15, rigorous=False)  
    
        P3=qs.sanunits._pumping.Pump('P103', ins=S1-0, outs='solvent103', P=101324, material='Stainless steel')
        P4=qs.sanunits._pumping.Pump('P104', ins=HX2-0, outs='additive103', P=101324, material='Stainless steel')
    
        DT1= Dissolution_tank('DT101', ins=(fs_stream.polymer, P3-0, P4-0), outs=('polymer_solution'), params=params, tau=24, kW_per_m3=14)
        P5=qs.sanunits._pumping.Pump('P105', ins=DT1-0, outs='polymer_ndegassolution', P=101324, material='Stainless steel')
        DG1 = Dissolution_degastank('DG101', ins=P5-0, outs=('polymer_degassolution'), tau=12, kW_per_m3=0.07, P=700)   

        #Bore solution
        M3= Mixer_water('M103', ins=fs_stream.water_boresol, outs='water_c', params=params, rigorous=False)
        HX3=qs.sanunits.HXutility('HX103', ins=M3-0, outs='water_h', T= params.T_boresol+273.15, rigorous=False)
        P6=qs.sanunits._pumping.Pump('P106', ins=HX3-0, outs='water101', P=101324, material='Stainless steel')
    
        P7=qs.sanunits._pumping.Pump('P107', ins=S1-1, outs='solvent104', P=101324, material='Stainless steel')
        DT2 = qs.sanunits.MixTank('DT102', ins=( P6-0, P7-0), outs=('bore_solution'), tau=1, kW_per_m3=1.25)  #Ranges from 1.0-1.5, which correspond to liquid-liquid mixing 
        P8=qs.sanunits._pumping.Pump('P108', ins=DT2-0, outs='bore_ndegassolution', P=101324, material='Stainless steel')
        DG2 = Bore_degasTank ('DG102', ins=(P8-0, fs_stream.nitrogen), outs=('bore_degassolution', 'nitrogen_out'), nitrogen_demand= 0.00056, params=params, tau=1, kW_per_m3=0.07)
        
        #Extrusion, coagulation & rinsing
        EX1= Extruder('EX201', ins=(DG1-0, DG2-0), outs=('lumen'), params=params,power_demand = 300)

        M4= Mixer_non_solvent('M201', ins=fs_stream.non_solvent, outs='non_solvent201', nonsolvent_per_membranearea=150, params=params, T=25+273, rigorous=False)
        HX4=qs.sanunits.HXutility('HX201', ins=M4-0, outs='nonsolvent_h', T= 20+273.15, rigorous=False)
        CB1= Coagulation_bath('CB201', ins=(EX1-0, HX4-0), outs=('precipitated_lumen', 'wastewater1'), solvent_out=0.7)
        
        if polymer_option in ('CA', 'CA_bioAA'):  
            RGT1= Regeneration_tank('RGT201', ins=(fs_stream.sodium_hydroxide, fs_stream.ethanol_rsolution, fs_stream.water_rsolution,'ethanolsol_recycled'), outs='regeneration_solution',
                 params=params, naoh_concentration=0.2, tau=0.5, kW_per_m3=2)   
            #HXB4= qs.sanunits.HXutility('HX201b', ins=RGT1-0, outs='regeneration_solutionh', T= 70+273.15, rigorous=False) #seems that regeneration is done at ambient temperature
            RB1= Regeneration_bath('RB201', ins=(CB1-0, RGT1-0), outs=('lumen_for_rinsing', 'wastewater1b'), tau=0.5, kW_per_m3=1.25) #from Ref 1 in regeneration tank
            SP1 = qs.sanunits.Splitter ('SP201', ins=RB1-1, outs=('purge', 3-RGT1), split=0.025)  #use a purge so Na-aacetate does not build up in the stream that is concentrated with EtOH

            stream_to_rinse = RB1-0
        else:
            stream_to_rinse = CB1-0
       
        M5= Mixer_rinsing_water('M202', ins=fs_stream.rinsing_water, outs='rinsing_water201', rinsing_per_membranearea=150, params=params, rigorous=False)
        HX5=qs.sanunits.HXutility('HX202', ins=M5-0, outs='rinsing_waterh', T= 50+273.15, rigorous=False)
        RT1= Rinsing_tank('RT201', ins=(stream_to_rinse, HX5-0), outs= ('lumen_after_rinsing', 'wastewater2'), solvent_out=0.99999, tau=0.5)
        
        #Preconditioning and drying
        M6= Mixer_glycerol_conditioning('M301', ins=fs_stream.glycerol_conditioning, outs='glycerol_conditioning301', fraction_filled=1.1, 
                                        params=params, rigorous=False)
        ST3=qs.sanunits.StorageTank('ST301', ins= M6-0, outs='glycerol_conditioning302', tau=336, 
                                    vessel_material='Stainless steel', vessel_type='Cone roof')  #low volatility
        P9=qs.sanunits._pumping.Pump('P301', ins=ST3-0, outs='glycerol_conditioning303', P=101324, material='Stainless steel')  #to the tank

        M7= Mixer_water_conditioning('M302', ins=fs_stream.water_conditioning, outs='water_conditioning301', 
                                    fraction_filled=1.1, params=params, rigorous=False)
        P10=qs.sanunits._pumping.Pump('P302', ins=M7-0, outs='water_conditioning301', P=101324, material='Stainless steel')  #to the tank

        DT3 = qs.sanunits.MixTank('MT301', ins=( P9-0, P10-0), outs=('conditioning_solutionc'), tau=0.5, kW_per_m3=1.25)  #Ranges from 1.0-1.5, which correspond to liquid-liquid mixing ##Residence time suggested for chta gpt for easy mixtures
        HX6=qs.sanunits.HXutility('HX301', ins=DT3-0, outs='conditioning_solutionh', T= 30+273, rigorous=False)
        P11=qs.sanunits._pumping.Pump('P303', ins=HX6-0, outs='conditioning_solution301', P=101324, material='Stainless steel')

        CT1 = Conditioning_tank('CT301', ins=(RT1-0, P11-0), outs=('conditioned_fiber'),tau=12, kW_per_m3=1.25, conditioning_temperature=30+273.5)   #high uncertainty in the residence time
        
        DY1 = Dryer('DY301', ins=(CT1-0), outs=('solid_fiber', 'wastewater3a'), final_moisture=0.0, drying_T = 50+273.15, dryer_efficiency = 0.7)
        HX7 = qs.sanunits.HXutility('HX302', ins=DY1-1, outs='wastewater3', T= 25+273.15, rigorous=False)
        
        #If we simulate the solids of PVP and PEG (we need to indicate a portion of this goes into wastewater streams!!, increased viscosity-fouling risk! not sure about precipitation)
        
        #Module assembly
        MD1 = Module_assembly('MD401', ins=(DY1-0, fs_stream.polysulfone_module, fs_stream.epoxy_module), outs=('module'), electricity_per_membranearea= 0.0667, 
                            epoxy_per_membranearea=0.127, params=params)
        
        #Waste water recycling and treatment
        if polymer_option in ('CA', 'CA_bioAA'):
            DT4 = qs.sanunits.MixTank('DT501', ins=( CB1-1, SP1-1, RT1-1, HX7-0), outs=('wastewater501'), tau=0.5, kW_per_m3=1.25)  #Ranges from 1.0-1.5, which correspond to liquid-liquid mixing 
        else:
            DT4 = qs.sanunits.MixTank('DT501', ins=( CB1-1, RT1-1, HX7-0), outs=('wastewater501'), tau=0.5, kW_per_m3=1.25)

        D1=qs.sanunits.ShortcutColumn('D501', ins=DT4-0, outs=('distillate', 1-M1), LHK=('H2O', 'NMP'), 
                                y_top=0.999, x_bot=0.001, k=3.05, is_divided=True, P=101324)   
        WT1 = Dilution_tank('WT501', ins=(D1-0, fs_stream.water_wwt), outs= fs_stream.diluted_wastewater, COD_target=425.34, tau=0.5, kW_per_m3=0.0)
        
        if polymer_option in ('CA', 'CA_bioAA'):
            sys_regen = qs.System('sys_regen', path=(RGT1, RB1, SP1), recycle=SP1-1)   # recycle back to RGT1

            sys_membrane=qs.System('sys_membrane', path=(M1, ST1, P1, HX1, S1, M2, ST2, P2, HX2, P3, P4, DT1, P5, DG1, M3, HX3, P6, P7, DT2, P8, DG2, EX1,M4, HX4, CB1, 
                                                         sys_regen, M5, HX5, RT1, M6, M7, P9, ST3, P10, DT3, HX6, P11, CT1, DY1, HX7, MD1, DT4, D1, WT1),recycle=D1-1)
      
        else:
            sys_membrane=qs.System('sys_membrane', path=(M1, ST1, P1, HX1, S1, M2, ST2, P2, HX2, P3, P4, DT1, P5, DG1, M3, HX3, P6, P7, DT2, P8, DG2, EX1,M4, HX4, CB1, 
                                                         M5, HX5, RT1, M6, M7, P9, ST3, P10, DT3, HX6, P11, CT1, DY1, HX7, MD1, DT4, D1,WT1), recycle=D1-1)
           
    else:
        raise RuntimeError('in function "system_membrane" argument "solvent_recycling" must be either "yes" or "no"')
    
    # RECYCLING NEEDS TO SOLVE FOR CA EVENTHOUGH MAYBE NO NEED RIGHT NOW AS WITH PSF WE CAN WORK IT OUT
    # Only recycling when water is less amount!

    sys_membrane.simulate()
    sys_membrane.show()
    sys_membrane.diagram()
    
    # ===========================================================================
    # Life cycle assessment (LCA)
    # ===========================================================================
   
    # --- Determine polymer impact item based on polymer_option ---

    if polymer_option == 'PSF':
        polymer_item = qs.ImpactItem.get_item('PSF_item')
    elif polymer_option == 'CA':
        polymer_item = qs.ImpactItem.get_item('CA_item')
    elif polymer_option == 'CA_bioAA':
        polymer_item = qs.ImpactItem.get_item('CA_bioAA_item')
    
    # --- Create the LCA ---
  
    lca = qs.LCA (system=sys_membrane, 
                lifetime=1, 
                lifetime_unit='year', 
                uptime_ratio=330/365,
                # Electricity: system power utility (kW) + all cooling HX units via electric chiller (COP=4)
                # and divide by 3600 kJ/hr -> kW; multiply by operating hours/year 
                # Steam: all heating HX units (net_duty > 0); multiply by operating hours/year 
                #(hasattr prevent a crash when it encounters units without it) 
                electricity_item = lambda: (sys_membrane.power_utility.rate
                                            + sum(abs(u.net_duty)
                                                    for u in sys_membrane.units
                                                    if hasattr(u, 'net_duty') and u.net_duty < 0
                                                    ) / 4 / 3600
                                            ) * (330 * 24),
                steam_item = lambda: (
                                        sum(
                                            u.net_duty
                                            for u in sys_membrane.units
                                            if hasattr(u, 'net_duty') and u.net_duty > 0
                                            ) / 1000 * (330 * 24)
                                    ),
                # --- Streams not auto-tracked by qsdsan ---
                PSF_module_item = lambda: fs_stream.polysulfone_module.F_mass * 24 * 330,
                epoxy_item = lambda: fs_stream.epoxy_module.F_mass * 24 * 330,
                water_wwt_item = lambda: fs_stream.water_wwt.F_mass * 24 * 330,
                NaOH_item = lambda: fs_stream.sodium_hydroxide.F_mass * 24 * 330,
                ethanol_item = lambda: fs_stream.ethanol_rsolution.F_mass * 24 * 330,
                water_rsolution_item = lambda: fs_stream.water_rsolution.F_mass * 24 * 330,
                )
        
    # Add untracked streams AFTER lca is created using add_other_item()
    lca.add_other_item(
        item       = polymer_item,
        f_quantity = lambda: fs_stream.polymer.F_mass * 24 * 330,
        unit       = 'kg'
        )
    # --- Add solvent only if recycling=yes (not auto-tracked in that case) ---
    if solvent_recycling == 'yes':
        lca.add_other_item(
            item       = qs.ImpactItem.get_item('NMP_item'),
            f_quantity = lambda: fs_stream.solvent.F_mass * 24 * 330,
            unit       = 'kg'
        )

    #--- Verification: Electricity and steam components --- (if needed, uncomment)
    # print('=== UNITS WITH net_duty > 0 (these count as STEAM, expec HX 101,102,103,202,301 and DY301) ===')
    # for u in sys_membrane.units:
    #        if hasattr(u, 'net_duty') and u.net_duty > 0:
    #         print(f'{u.ID}: net_duty = {u.net_duty:.2f} kJ/hr')

    # print('\n=== UNITS WITH net_duty < 0 (cooling — HX201 and HX302 expected here) ===')
    # for u in sys_membrane.units:
    #     if hasattr(u, 'net_duty') and u.net_duty < 0:
    #         print(f'{u.ID}: net_duty = {u.net_duty:.2f} kJ/hr')    
      
    # lca.show()

    # --- Parameter needed to estimate LCA results per FU (m^2) ---
    polymer_flow = fs_stream.polymer.F_mass  #kg/h
    membrane_area_per_hour = polymer_flow / params.polymermass_to_membranearea   #m^2/h
    membrane_area_per_year = membrane_area_per_hour*24*330
    
    
    return sys_membrane, lca, membrane_area_per_year

system, lca, membrane_area_per_year= system_membrane(8000, 'PSF', 'no')
print(membrane_area_per_year)
#system.save_report(file='test_recycling_150.xlsx')

LCAresults = LCAResults(lca=lca, 
                        membrane_area_per_year=membrane_area_per_year, 
                        polymer_option='PSF', 
                        solvent_recycling='no'
                        )


LCAresults.show()
#LCAresults.save()

# ===========================================================================
# Uncertainty analysis
# ===========================================================================

# uncertainty, stats = uncertainty_analysis(system, 'all', 3000, membrane_area_per_year,
#                        'PSF', 'yes')

# # # Access full results table if needed
# #only main metics
# # --- Final metric scores ONLY (what you report in paper) ---
# actual_metric_names = [
#     f'{m.name} [{m.units}]'
#     for m in uncertainty.metrics
#     if m.element == 'LCA'
# ]
# final_scores = uncertainty.table.loc[:, 'LCA'][actual_metric_names]
# final_scores.to_csv('uncertainty_PSF_yes_3000.csv')

#lca bacground and results
#results_table = uncertainty.table.loc[:, 'LCA']
#whole table
#results_table.to_csv('uncertainty_PSF_no_1000.csv')

# ===========================================================================
# Sensitivity analysis
# ===========================================================================

# sensitivity_model, r_df, p_df = sensitivity_analysis(
#     sys                    = system,
#     num_samples            = 3000,
#     membrane_area_per_year = membrane_area_per_year,
#     polymer_option         = 'PSF',
#     solvent_recycling      = 'yes',
#     threshold              = 0.1,   # exclude parameters with |rho| < 0.1 everywhere
# )



# %%
# ===========================================================================
# How to get the results per m^2 membrane (results from LCA are per year)
# ===========================================================================
#Get LCA results
# impacts_per_m2 = {
#     indicator: value / membrane_area_per_year
#     for indicator, value in lca.get_total_impacts().items()
# }

# print('\n=== LCA RESULTS PER m2 MEMBRANE ===')
# for indicator, value in impacts_per_m2.items():
#     print(f'  {indicator}: {value:.6e}')

# # ================================================================
# # 2. CONTRIBUTION BY STREAM
# # ================================================================
# stream_contributions1 = lca.get_impact_table('Stream')
# stream_contributions2 = lca.get_impact_table ('Other')
# print(stream_contributions1)
# print(stream_contributions2)


# #See all streams the LCA is tracking with their annual flows (this can help us to identify if we do not need to add manually in the LCA)
# print('=== Streams tracked by LCA ===')
# for s in lca.system.streams:
#     if hasattr(s, 'stream_impact_item') and s.stream_impact_item is not None:
#         item = s.stream_impact_item
#         print(f'{s.ID:30s} | F_mass={s.F_mass:.4f} kg/hr | item={item.ID}')


# %%
