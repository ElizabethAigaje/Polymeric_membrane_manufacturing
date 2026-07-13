"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 12/05/2025
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom
import math


class Extruder(qs.SanUnit):
    """ Creates a extrusion/spinning system for fiber formation.

    Parameters
    ----------
    ID : str
        ID of extruder.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    params: str
        Variable that stores system-level parameters, [-] (see class 'MembraneParams' in the system-membrane-manufacutring.py file)
    power_demand: str
        Power demand of a extruder that produces 4000000 modules per year, [kW] (defualt value is 300 kW)


    Notes  
    -----
    The extrusion/spinning system includes the extruder and pumps for bore and polymer solutions, power demand includes extruder and roller (winding) systems [1].
    For the cost calculation, it is assumed that the selected design is the biggest one, and any increase in the capacity adds another extrusion system.
    We assumed the base capacity is 400000 modules annualy as an intermediate capacity. 

    References
    ----------
    [1] HFSM. Industrial scale hollow fiber spinning machine. https://hollowfiberspinningmachine.com/

    """

    #Information in the design results  ##CHECK THIS WITH THE DESIGN OF THE EXTRUDER
    _units={'Electricity_demand':'kwh-h', 'Capacity': 'modules/year'}      

    def __init__(self, ID, ins=(), outs=(), thermo=None, power_demand=300, params=None):
        super().__init__(ID, ins, outs,thermo)  
    
        self.params=params
        self.power_demand = power_demand


    _N_ins=2
    _N_outs=1

    def _run(self):

        polymer_solution, bore_solution = self.ins           
        lumen, = self.outs               

        lumen.imass['NMP']=polymer_solution.imass['NMP']+bore_solution.imass['NMP']
        lumen.imass['PEG']=polymer_solution.imass['PEG']
        lumen.imass['polysulfone']=polymer_solution.imass['polysulfone']
        lumen.imass['H2O']=bore_solution.imass['H2O']
 
        lumen.phase ='l'

    def _design(self):

        p=self.params
        base_extruder_capacity = 400000 
        base_power= self.power_demand   # 300 #kW  this power includes the pumping, spinning, coagulation and winding.
        
        polymer_needed=(p.v_polymersol * p.d_polymersol *
                         p.polymer_fraction) / 1000.0 / 24.0  #kg/h
        self.needed_capacity= (polymer_needed/p.polymermass_to_membranearea)*365*24/p.membranearea_per_module
        
         
        electricity_demand = (self.needed_capacity/base_extruder_capacity)*base_power  #kWh #gonna exclude the math ceil to just have proportional

        self.design_results['Electricity_demand'] = electricity_demand       
        self.design_results['Capacity'] = self.needed_capacity                           
        self.add_power_utility(electricity_demand)    #not used the 180 kW because it include coagulation, washing, drying and widing and we are using just spnining

    def _cost(self):
        
        base_extruder_capacity = 400000  # modules per year    
        base_purchase_cost= 1000000   # function of bundles capacity #Need this number.
                                #Needs more information https://hollowfiberspinningmachine.com/
                                # https://www.me-sep.com/me-sep-consultancy-research/
        purchase_cost = base_purchase_cost*math.ceil(self.needed_capacity/base_extruder_capacity)   #We will not apply economy of scale, since the n factor is not characterize for this type of equipment 
                                                                                    #and we will assume this is the biggest capacity
        self.baseline_purchase_costs['Extruder'] = purchase_cost   
        self.F_D['Extruder']=self.F_P['Extruder']=self.F_M['Extruder']=1  #design, pressure, material factors
