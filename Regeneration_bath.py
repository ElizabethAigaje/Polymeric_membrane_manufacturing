"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 03/02/2026
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Regeneration_bath(qs.sanunits.MixTank): 
    """ Creates a regeneration bath system with NaOH solution.

    Parameters
    ----------

    ID : str
        ID of the equipment.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    tau : float
        Regeneration time, [h] 
    kW_m3 : float
        Power required for mixing (the default is 1.25 for liquid-liquid mixing [1])

    Notes
    -----
    Stirred tank for the regeneration bath following the coagulation and rinsing tank designs.
    Residence time comes from industry input in continous processes.

    References
    ---------
    [1] Sinnott, R. K., & Towler, G. P. (2020). Chemical engineering design (6th edition). Elsevier, pp. 625.
        
    """

    _N_ins=2
    _N_outs=2


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                  tau=(), kW_per_m3=1.25):   #liquid-liquid mixing
        super().__init__(ID, ins, outs, Thermo) 
      
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):
            
        # Rinsing to remove the solvent remaining in the lumen
        lumen_precipitated, regeneration_solution = self.ins           
        fiber_out, rsolution_out= self.outs               

        fiber_out.imass['polysulfone']=lumen_precipitated.imass['polysulfone']
        fiber_out.imass['PEG']=lumen_precipitated.imass['PEG']
        fiber_out.imass['H2O']= lumen_precipitated.imass['H2O']
        fiber_out.imass['NMP']= lumen_precipitated.imass['NMP']

        ###
        rsolution_out.copy_like(regeneration_solution)
        rsolution_out.imass['H2O']= regeneration_solution.imass['H2O']
        rsolution_out.imass['NaOH']= regeneration_solution.imass['NaOH']


        fiber_out.phase ='l'    #even though fiber is in solid, water and solvent are in liquid state 
        rsolution_out.phase ='l'
            

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   