"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 12/01/2025
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Bore_degasTank(qs.sanunits.MixTank): 
    """ Creates an agitated tank for degassification of the bore solution with nitrogen.

    Parameters
    ----------
    ID : str
        ID of the dissolution tank.
    ins : Iterable(stream)
        polymer, solvent, additive.
    outs: Iterable(stream)
        dissolution.
    nitrogen_demand: 
        The nitrogen demand bor degassing bore solution, [L/m^2 membrane] (the default is 0.00056 [1])
    params:
        Variable that stores system-level parameters, [-] (see class 'MembraneParams' in the system-membrane-manufacutring.py file).
    tau : float
        Dissolution time, [hr]
    kW_m3 : float
        Power required for degassing (the default is 0.07 (ranges from 0.04-0.1), which correspond to mild blending, mixing [1])

    Notes
    -----
    The nitrogen mass flow is estimated from the polymer mass flow and polymer mass per m^2 of membrane.

    References
    ---------
    [1] Yadav, P., Ismail, N., Essalhi, M., Tysklind, M., Athanassiadis, D., & Tavajohi, N. (2021).
        Assessment of the environmental impact of polymeric membrane production. Journal of Membrane Science, 622, 118987. 
        https://doi.org/10.1016/j.memsci.2020.118987
    [2] Sinnott, R. K., & Towler, G. P. (2020). Chemical engineering design (6th edition). Elsevier, pp. 625.
    
    """

    _N_ins=2
    _N_outs=2   


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 nitrogen_demand= 0.00056, params=None, tau=(), kW_per_m3=0.07):   #max solubility at 120 C for LPDE in solvent, power required to mix slurry suspensions 
        super().__init__(ID, ins, outs, Thermo) 
        self.nitrogen_demand = nitrogen_demand
        self.params=params
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):

        p=self.params

        """Calculate nitrogen demand."""
        bore, nitrogen = self.ins           
        boredegas, nitrogen_out = self.outs               
        
        # Calculate the amount of polymer required 
        
        boredegas.copy_like(bore) 

        
        polymer_needed=(p.v_polymersol * p.d_polymersol *
                         p.polymer_fraction) / 1000.0 / 24.0  #kg/h 
        nitrogen_needed = (polymer_needed/p.polymermass_to_membranearea)*self.nitrogen_demand*1.132/1000     #nitrogen density at 25 C and 1 bar
        
        nitrogen.imass['N2']=nitrogen_needed
      
        
        nitrogen_out.imass['N2']=nitrogen_needed
      
 
        nitrogen.phase = nitrogen_out.phase ='g'
        boredegas.T= bore.T 
            

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   
