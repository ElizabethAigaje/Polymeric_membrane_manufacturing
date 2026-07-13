"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 2/4/2026
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom
from thermosteam import Chemical
import numpy as np

qs.CEPCI=567.5  #CEPCI 2017 IN QSDSAN

class Module_assembly(qs.SanUnit):
    """ Creates a module to estimate materials and energy flows during boundling (membrane formation) and module assembly  

    Parameters  #NEED TO BE DONE
    ----------
    ID : str
        ID of the unit
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    electricity_per_membranearea: float, optional
        Electricity demand to assembly the module per area of membrane [kWh/m^2] (the default is 0.0667 kWh/m^2 [1])
    epoxy_per_membranearea: float, optional
        Mass of epoxy resin needed for the cartridge per area of membrane [kg/m^2] (the default is 0.127 kg/m^2 [1])
    params:
        Variable that stores system-level parameters, [-] (see class 'MembraneParams' in the system-membrane-manufacutring.py file).

    Notes
    -----
    This module is based on energy and mass flows reported in reference [1]. For the hollow fiber wound, a mechanical arm is assumed to assembled the fibers [1].
    Energy demand for boundling, gluing and cutting is 6.7E-3, 3.6E-2 and 2.4 E-2 kWh per m^2, respectively, with a total of 0.0667 kWh per m^2 membrane.
    For the module materials, we consider:  the housing, end caps and permeate ports to be made from polysulfone as indicated in the cytiva manual. The amount 
    required will be approximated from the cartridge dimensions assuming the thickness (as it is not given). The amount of the epoxy fiber potting
    will be taken from the epoxy resin inventory from reference [1].
    
    !!!COST OF THE ARM ??
    References
    ----------
    [1] F. Prezelus, et al. A generic process modelling-LCA approach for UF membrane fabrication: Application to cellulose acetate membranes. 
    Journal of Membrane Science 618 (2021) 118594.

    """

    _N_ins=3
    _N_outs=1   

    _units={'Electricity_demand':'kW'}   #units for the design results

    def __init__(self, ID, ins=(), outs=(), thermo=None, electricity_per_membranearea= 0.0667, 
                 epoxy_per_membranearea=0.127,
                 params=None ):


        super().__init__(ID, ins, outs, thermo)    #inherinting the funtionalities of Units
        self.electricity_per_membranearea=electricity_per_membranearea
        self.epoxy_per_membranearea=epoxy_per_membranearea
        self.params=params


    def _run(self):
        p=self.params

        """ Calculate the amount f PSF for the module and its components and the amount of adhesive. """
        feed, psf, epoxy = self.ins     
        module, = self.outs

    # Copy everything to dry stream first
        module.copy_like(feed)

        self.polymer_needed=(p.v_polymersol * p.d_polymersol *
                         p.polymer_fraction) / 1000.0 / 24.0  #kg/h of polymer
        psf_needed = (self.polymer_needed/p.polymermass_to_membranearea)*p.cartridge_mass_per_area    #kg/h of psf for cartridge

        epoxy_needed= (self.polymer_needed/p.polymermass_to_membranearea)*self.epoxy_per_membranearea 

        psf.empty()
        psf.imass['polysulfone'] = psf_needed

        epoxy.empty()
        epoxy.imass['epoxy'] = epoxy_needed  

        module.copy_like(feed)   #outlet will have everyhting that comes in plus psf and epoxy!
        module.imass['polysulfone']=feed.imass['polysulfone']+psf.imass['polysulfone']
        module.imass['epoxy']=epoxy.imass['epoxy']

        module.T=25+273      #Room temperature for the module
        module.phase ='s'      #Lets see how this work out

       
    def _design(self):

        p=self.params
        """ Calculate the electricity demand and store in design results. """

        electricity_demand= (self.polymer_needed/p.polymermass_to_membranearea)*self.electricity_per_membranearea #kWh/h
        self.design_results['Electricity_demand']=electricity_demand
        self.add_power_utility(electricity_demand) 

    ## I may need to work on this cost for the electric arm
    # def _cost(self):
    #     """ Calculate purchase cost based on evaporation rate. """

    #     S=self.design_results['Evaporation_rate']
    #      # Avoid log(0)
    #     S = max(S, 1e-6)

    #     if S<=3000:      #maximum size of a spray dryer taken as reference  
    #         lnS= np.log(S)
    #         purchase_cost=np.exp(8.5133 + 0.9847*lnS-0.0561*(lnS**2))*(qs.CEPCI/567)   
    #     else:
    #         number_dryers= np.ceil(S/3000)
    #         lnS= np.log(3000)
    #         purchase_cost_unit=np.exp(8.5133 + 0.9847*lnS-0.0561*(lnS**2))*(qs.CEPCI/567) 
    #         purchase_cost=number_dryers*(purchase_cost_unit)
        
    #     self.baseline_purchase_costs['Dryer'] = purchase_cost  
    #     self.F_D['Dryer']=self.F_P['Dryer']=self.F_M['Dryer']=1  
    #     self.F_BM['Grinder']=2.06

