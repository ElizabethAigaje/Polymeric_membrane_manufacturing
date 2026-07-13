"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 1/27/2026
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Conditioning_tank(qs.sanunits.MixTank): 
    """ Creates a conditioning bath system with a glycerol solution prior drying.

    Parameters
    ----------

    ID : str
        ID of conditioning bath.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    tau : float
        Rinsing time, [h] 
    kW_m3 : float
        Power required for mixing (the default is 1.25 for liquid-liquid mixing [1])
    conditioning_temperature :
        Conditioning temperature, [K] (the default is 30 C [2])

    Notes
    -----
    This step is required to alliviate membrane shrinkage and prevent collapsing during drying [3]. Soaking the fibers into a conditioning bath is required. 
    We assume the glycerol solution stays in the fiber [2] (and the glycerol remains until it is first used [4]) and therefore upstream, the demand of glycerol solution
    is calculated based on the membrane pore volume. Similarly to raising tanks, we assume a continous (agitated) tank that is already filled with the solution 
    (we assume there is no contamination [3]). This model acccounts for the preconditioning materials demand. Moreover, the water from the fiber before conditioning in the pores 
    is assumed to stay in the membrane for simplicity in the analysis. 
    
    References
    ---------
   
    [1] Sinnott, R. K., & Towler, G. P. (2020). Chemical engineering design (6th edition). Elsevier, pp. 625.
    [2] Prezelus, F.; Tiuta-Barna, L.; Guigui, C.; and Remigy, J. A generic process modelling-LCA approach for UF membrane fabrication:
    Application to cellulose acetate membranes. Journal of Membrane Science 618 (2021) 118594
    [3] Gao, J., & Chung, T. Influence of contaminants in glycerol/water mixtures during post-treatment  on 
    physicochemical properties and separation performance of air-dried  membranes. Journal of Membrane Science 572 (2019) 223–229
    [4] CYTIVA. Operating handbook: Hollow fiber cartridges for membrane separations, p.9.   
    """

    _N_ins=2
    _N_outs=1


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 tau=(), kW_per_m3=1.25, conditioning_temperature = 30+273.15):   
        super().__init__(ID, ins, outs, Thermo) 
      
        self.tau=tau
        self.kW_per_m3=kW_per_m3
        self.conditioning_temperature = conditioning_temperature

    def _run(self):
            
        # Rinsing to remove the solvent remaining in the lumen
        fiber_in, preconditioning_in = self.ins           
        fiber_conditioned_out, = self.outs               

        fiber_conditioned_out.imass['polysulfone']=fiber_in.imass['polysulfone']
        fiber_conditioned_out.imass['PEG']=fiber_in.imass['PEG']
        fiber_conditioned_out.imass['H2O']= fiber_in.imass['H2O'] + preconditioning_in.imass['H2O']
        fiber_conditioned_out.imass['NMP']= fiber_in.imass['NMP']
        fiber_conditioned_out.imass['glycerol']= preconditioning_in.imass['glycerol']
        
 
        fiber_conditioned_out.phase ='l'    #even though fiber is in solid, water and solvent are in liquid state 
        fiber_conditioned_out.T= self.conditioning_temperature
            

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   