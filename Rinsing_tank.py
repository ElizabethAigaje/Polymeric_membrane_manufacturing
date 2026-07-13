"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 12/9/2025
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Rinsing_tank(qs.sanunits.MixTank): 
    """ Creates a rinsing bath system with water.

    Parameters
    ----------

    ID : str
        ID of extruder.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    solvent_out :
        Mass fraction of solvent difussing out in the rinsign stage, (the default is 0.9999)
    tau : float
        Rinsing time, [h] 
    kW_m3 : float
        Power required for mixing (the default is 1.25 for liquid-liquid mixing [1])

    Notes
    -----
    Stirred tank for the rinsing bath following description in reference [2].
    The rinsing bath requires longer residence time or higher bath temperature. We will assume it is performed at 50 C [3] and use the residence 
    time suggested by the industry members.
    The costs will be calculated separetely from the extruder and coagulation bath system, as rinsing conditions are specific to this system.

    References
    ---------
    [1] Sinnott, R. K., & Towler, G. P. (2020). Chemical engineering design (6th edition). Elsevier, pp. 625.
    [2] Prezelus, F.; Tiuta-Barna, L.; Guigui, C.; and Remigy, J. A generic process modelling-LCA approach for UF membrane fabrication:
    Application to cellulose acetate membranes. Journal of Membrane Science 618 (2021) 118594
    [3] Chidambaran, R.; Shara, D.; Raina, P.; Das, S. Preapration of high performance ultra filtation hollow fiber membrane. US Patent
    US2012/0111790 A1 (2012)
        
    """

    _N_ins=2
    _N_outs=2


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 solvent_out=0.9999, tau=(), kW_per_m3=1.25):   #max solubility at 120 C for LPDE in solvent, power required to mix slurry suspensions 
        super().__init__(ID, ins, outs, Thermo) 
      
        self.solvent_out = solvent_out
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):
            
        # Rinsing to remove the solvent remaining in the lumen
        lumen_precipitated, rinsing_solution = self.ins           
        fiber_out, rinsing_out= self.outs               

        fiber_out.imass['polysulfone']=lumen_precipitated.imass['polysulfone']
        fiber_out.imass['PEG']=lumen_precipitated.imass['PEG']
        fiber_out.imass['H2O']= lumen_precipitated.imass['H2O']
        fiber_out.imass['NMP']= (1-self.solvent_out)*lumen_precipitated.imass['NMP']

        rinsing_out.imass['H2O']= rinsing_solution.imass['H2O']
        rinsing_out.imass['NMP']=self.solvent_out*lumen_precipitated.imass['NMP']
 
        fiber_out.phase ='l'    #even though fiber is in solid, water and solvent are in liquid state 
        rinsing_out.phase ='l'
            

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   