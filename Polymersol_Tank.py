"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 11/27/2025
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Dissolution_tank(qs.sanunits.MixTank): 
    """ Creates an agitated tank to prepare the polymer solution.

    Parameters
    ----------
    ID : str
        ID of the dissolution tank.
    ins : Iterable(stream)
        polymer, solvent, additive.
    outs: Iterable(stream)
        dissolution.
    params:
        Variable that stores system-level parameters, [-] (see class 'MembraneParams' in the system-membrane-manufacutring.py file).
    tau : float
        Dissolution time, [hr]
    kW_m3 : float
        Power required for mixing (the default is 14 [1])

    Notes
    -----
    The total polymer requirement is calculated from the `v_polymersol`, 'polymer_fraction' and the `rho_polymersol` stored in the params variable.
    For the energy required we utilize literature data for dissoling polymers into solvents. Ref [1] reports that 0.8 kWh is required to
    dissolve 1 lb of polyolefins. This translate to 14 kwh per m^3. Reference [2] calculates approximately 22 kwh per m^3 (for this last, the tank volume
     is not reported, therefore this estimation is based on our tank volume). We selected data from ref 1 as our baseline.

    References
    ---------
    [1] Naviroj, P., Treacy, J. and Urffer, C. Chemical recycling of plastics by dissolution. University of Pennsylvania (2019)
    [2] Prezelus, F., et al. A generic process modelling – LCA approach for UF membrane fabrication: Application to cellulose acetate membranes
        Journal of Membrane Science 618 (2021) 118594.
        
    """

    _N_ins=3
    _N_outs=1


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 params=None, tau=(), kW_per_m3=14):   #max solubility at 120 C for LPDE in solvent, power required to mix slurry suspensions 
        super().__init__(ID, ins, outs, Thermo) 
        self.params = params
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):

        p=self.params

        """Calculate material flow distribution and phase conditions."""
        polymer, solvent, additive = self.ins           
        dissolution, = self.outs               
        
        # Calculate the amount of polymer required 

        polymer_needed=(p.v_polymersol * p.d_polymersol *
                         p.polymer_fraction) / 1000.0 / 24.0  #kg/h 
        polymer.imass['polysulfone']=polymer_needed
       
        solvent_mass=solvent.imass['NMP']
        additive_mass=additive.imass['PEG']
        
        dissolution.imass['NMP']=solvent_mass
        dissolution.imass['PEG']=additive_mass
        dissolution.imass['polysulfone']= polymer_needed
 
        dissolution.phase ='l'
        dissolution.T= solvent.T 
            

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   