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
import math
from biosteam.units.design_tools.specification_factors import vessel_material_factors
from biosteam.units.design_tools.tank_design import (
    TankPurchaseCostAlgorithm,
    compute_number_of_tanks_and_purchase_cost,
    storage_tank_purchase_cost_algorithms,
    mix_tank_purchase_cost_algorithms)


class Dissolution_degastank(qs.sanunits.MixTank): 
    """ Creates an agitated tank under vacuum for degasification of the polymer solution.

    Parameters
    ----------
    ID : str
        ID of the dissolution tank.
    ins : Iterable(stream)
        polymer, solvent, additive.
    outs: Iterable(stream)
        dissolution.
    P : str
        Absolute operating vacuum pressure assuming a barometric pressure of 760 torr, [torr] (the default is 700 torr for rough vacuum region [1], ranges from 700 to 760)
    tau : float
        Dissolution time, [hr]
    kW_m3 : float
        Power required for degassing (the default is 0.07 (ranges from 0.04-0.1), which correspond to mild blending, mixing [2])

    Notes
    -----
    This tank includes mild agitation under a vacuum environemnt. The cost of the vacuum system for degasing is included in the tank cost.
    The cost logarithm is for a liquid-ring pump, taken from reference [3]

    References
    ---------

    [1] Chidambaran, R., Sharma, D., Raina, P., Das, S. Preparation of high performance ultra filtration hollow fiber membrane.
        Patent US 2012/0111790 A1
    [2] Sinnott, R. K., & Towler, G. P. (2020). Chemical engineering design (6th edition). Elsevier, pp. 625.
    [3] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Purchase Costs of Other Chemical Processing Equipment. In Product
        and Process Design Principles; Wiley, 2017; pp 479-485.

    
    """

    _N_ins=1
    _N_outs=1


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 P = 700,  tau=(), kW_per_m3=0.07):   
        super().__init__(ID, ins, outs, Thermo) 
        self.P=P
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):

        self.outs[0].copy_like(self.ins[0])  
            

    # def _design(self):
    #     """Inherits design calculation from MixTank."""
    #     qs.sanunits.MixTank._design(self)
    #     #Accounts for the number of tanks needed to work in paralell to cover continous process

    # def _cost(self):
    #     """Inherits cost calculation from MixTank."""
    #     qs.sanunits.MixTank._cost(self)
       
         
    @property
    def vessel_type(self):
        return self._vessel_type
    @vessel_type.setter
    def vessel_type(self, vessel_type):
        if vessel_type in self.purchase_cost_algorithms:
            self._vessel_type = vessel_type
            self.purchase_cost_algorithm = self.purchase_cost_algorithms[vessel_type]
        else:
            raise ValueError(f"vessel type '{vessel_type}'"
                             "is not avaiable; only the following vessel "
                             "types are available: "
                            f"{', '.join(self.purchase_cost_algorithms)}")
    
    def _design(self):
        design_results = self.design_results
        design_results['Residence time'] = tau = self.tau
        design_results['Total volume'] = tau * self.F_vol_out / self.V_wf


    def _cost(self):
        design_results = self.design_results

        V = design_results['Total volume']
        W = 5 + (0.0298+0.03088*(math.log(self.P))-0.0005733*((math.log(self.P)))**2)*((V*35.3147)**0.66)   #self.P in torr and V from m^3 in ft^3
        cost_pump= (8250*(W)**0.37)*567.5/567   #  CEPCI updated since biosteam costs are in 2017 and Seider book is CE=567
        
        N, Cp = compute_number_of_tanks_and_purchase_cost(
            V, self.purchase_cost_algorithm
        )
        if N:
            self.parallel['self'] = N
            default_material = self.purchase_cost_algorithm.material
            self.baseline_purchase_costs['Tank']  = (Cp / vessel_material_factors.get(default_material, 1.))+N*cost_pump
            self.add_power_utility(self.kW_per_m3 * V / N)
        


   