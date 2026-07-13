"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 03/25/2026
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Regeneration_tank(qs.sanunits.MixTank): 
    """ Creates a agitated tank to prepare the regeneration solution.

    Parameters
    ----------
    ID : str
        ID of the dissolution tank.
    ins : Iterable(stream)
        sodiumhydroxide, water.
    outs: Iterable(stream)
        dissolution.
    params:
        Variable that stores system-level parameters, [-] (see class 'MembraneParams' in the system-membrane-manufacutring.py file).
    naoh_concentration:
        NaOH concentration, [M] (the default is 0.2 M) [1]
    tau : float
        Dissolution time, [hr]
    kW_m3 : float
        Power required for mixing (the default for suspension of solid particles is 2 [2])

    Notes
    -----
    The total polymer requirement is calculated from the `v_polymersol`, 'polymer_fraction' and the `rho_polymersol` stored in the params variable.
    This is the dissolution tank to prepare the regeneration solution, including the recycling stream.

    References
    ---------
    [1] He, X. Optimization o deacetylation process for regenerated cellulose hollow fiber membranes. 
        International Journal of Polymer Science Volume 2017, Article ID 3125413, 8 pages.
    [2] Seider, W. D.; Seader, J. D.; Lewin, D. R. PRODUCT & PROCESS DESIGN PRINCIPLES: SYNTHESIS, ANALYSIS AND EVALUATION, 
        (With CD); John Wiley & Sons, 2009, p. 470
        
    """
    #CHANGE ALL THIS

    _N_ins=4
    _N_outs=1


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 params=None, 
                 naoh_concentration=0.2, tau=(), kW_per_m3=2):   #max solubility at 120 C for LPDE in solvent, power required to mix slurry suspensions 
        super().__init__(ID, ins, outs, Thermo) 
        self.params=params
        self.naoh_concentration=naoh_concentration
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):
        p=self.params

        """Calculate material flow distribution and phase conditions."""
        sodiumhydroxide, ethanol, water, recycle_stream = self.ins           
        dissolution, = self.outs               
        
        # Calculate inlet flows 
        
        polymer_needed=(p.v_polymersol * p.d_polymersol *
                        p.polymer_fraction) / 1000.0 / 24.0  #kg/h 

        #Solution volume needed
        AGU_MW=162   #NaOH -->Na-acetate   kg/kmol
        kmol_AGU= polymer_needed/AGU_MW    #kmol/h moles of anhydride glusoe unit in the celulose acetate
        CA_DS = 2.25  #Range of DS is about 2-2.5

        kmol_acetylgroups= kmol_AGU*CA_DS        #kmol of acetyl groups due to acetylation in cellulose acetate 
        #mol_acetlygrous = mol_NaOH  for regneration

        naoh_needed = kmol_acetylgroups*40     #kg/h NaOH are consumed during regeneration

        solution_volume = (naoh_needed/40)*1000/self.naoh_concentration   # solution volume L
        solution_needed = solution_volume * 789/1000    #kg/h, assuming most of the solution is ethanol 789 kg/m^3

        ethanol_needed= (solution_needed - naoh_needed)*0.96   #kg EtOH/h  #96% aquous ethanol ref XX
        water_needed= (solution_needed - naoh_needed)*0.04   #kg EtOH/h

        #Then inlets are:

        sodiumhydroxide.imass['NaOH']=naoh_needed
        ethanol.imass['EtOH']=ethanol_needed-recycle_stream.imass['EtOH']
        water.imass['H2O'] = water_needed-recycle_stream.imass['H2O']

        #Outlets are: 
       
        naoh_mass = sodiumhydroxide.imass['NaOH']
        ethanol_mass= ethanol.imass['EtOH'] + recycle_stream.imass['EtOH']
        water_mass = water.imass['H2O'] + recycle_stream.imass['H2O']
        
        dissolution.imass['NaOH']= naoh_mass
        dissolution.imass['EtOH']= ethanol_mass
        dissolution.imass['H2O']= water_mass
 
        dissolution.phase ='l'
        dissolution.T= 25+273
            
        # print(self.ID, "RUNNING")
        # print("NaOH kg/h =", naoh_mass)
        # print("EtOH kg/h =", ethanol_mass)
        # print("Water kg/h =", water_mass)

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   