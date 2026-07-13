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
from thermosteam import Chemical
import numpy as np

qs.CEPCI=567.5  #CEPCI 2017 IN QSDSAN

class Dryer(qs.SanUnit):
    """ Creates a dryer unit to remove water from membrane fibers

    Parameters:
    ----------
    ID : str
        ID of the dryer
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    final_moisture : float, optional
        Moisture content of the membrane fibers, (the default is 0).
    drying_T : float, optional
        Drying temperature, [K] (the default is 323 K [1]).
    drying_efficiency : float, optional
        Drying efficiency, (the default is 0.7).

    Notes
    -----
    The dryer is modeled as a continuous convective oven, consistent with industrial hollow-fiber manufacturing lines. 
    Drying was modeled as steady-state removal of residual water to a target moisture content, with thermal energy demand calculated from 
    sensible heating and latent heat of vaporization. Cost equation is based on a spray dryer to assume continous operation
    (assuming a hot-air tunnel/conveyor dryer, which are calculated based on the evaporation rate) [2].

    References
    ----------
    [1] Patent WO 2004/089520 A1 (2003) 
    [2] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Purchase Costs of Other Chemical Processing Equipment. In Product
        and Process Design Principles; Wiley, 2017; pp 484.

    """

    _N_ins=1
    _N_outs=2   # Dry fiber and evaporated water

    _units={'Evaporation_rate':'lb/h', 'Duty':'kJ/h'}   #units for the design results

    def __init__(self, ID, ins=(), outs=(), final_moisture=0., drying_T = 50+273.15, dryer_efficiency = 0.7, thermo=None):
        super().__init__(ID, ins, outs, thermo)    #inherinting the funtionalities of Units
        self.final_moisture=final_moisture 
        self.drying_T = drying_T
        self.dryer_efficiency = dryer_efficiency

    def _run(self):
        """ Calculate the amount of water out in the product and the vapor stream. """
        feed, = self.ins
        dry_fiber, vapor = self.outs
    # Copy everything to dry stream first
        dry_fiber.copy_like(feed)
        vapor.empty()

        # Masses
        m_water = feed.imass['H2O']
        m_total = feed.F_mass
        m_solids = m_total - m_water     #Fiber, solvent,glycerol, glycol

        # Target outlet water (wet basis)
        m_water_out = self.final_moisture/(1 - self.final_moisture) * m_solids #water content in the final solid fiber

        # Evaporated water
        m_evap = max(m_water - m_water_out, 0)

        # Update outlet streams
        dry_fiber.imass['H2O'] = m_water_out
        vapor.imass['H2O'] = m_evap

        # Everything else stays with dry fibers
        for chem in feed.chemicals:
            if chem.ID != 'H2O':
                vapor.imass[chem.ID] = 0

        # Set temperatures
        dry_fiber.T = self.drying_T
        vapor.T = self.drying_T
        vapor.phase = 'l'

        # Save for design step
        self.m_evap = m_evap   #kg/h
        
    def _design(self):
        """ Calculate evaporation rate, utility demand and store in design results. """

        feed = self.ins[0]
        evap_rate = self.m_evap*2.2    #lb/h for costing
        self.design_results['Evaporation_rate']=evap_rate

        #--- Energy calculation ---

        # Sensible heat 
        cp= 3e3    #J/kg/k needs an approximation for fibre wet LOOK AT THIS!
        dT = self.drying_T - feed.T
        Q_sens = feed.F_mass*cp*dT #J/h

        # Latent heat
        Water = Chemical('Water', cache=True)
        H_vap = Water.Hvap(self.drying_T) #J/mol
        Q_lat = self.m_evap*H_vap*1000/18 #J/h

        Q_total = Q_sens + Q_lat   #J/h

        #Accounting for dryier efficiency

        Q_required = (Q_total/self.dryer_efficiency)/1000  #kJ/h
    
        self.design_results['Duty']=Q_required                              
        self.add_heat_utility(Q_required, self.drying_T)    #do not add the agent so the system defines the more suitable

    def _cost(self):
        """ Calculate purchase cost based on evaporation rate. """

        S=self.design_results['Evaporation_rate']
         # Avoid log(0)
        S = max(S, 1e-6)

        if S<=3000:      #maximum size of a spray dryer taken as reference  
            lnS= np.log(S)
            purchase_cost=np.exp(8.5133 + 0.9847*lnS-0.0561*(lnS**2))*(qs.CEPCI/567)   
        else:
            number_dryers= np.ceil(S/3000)
            lnS= np.log(3000)
            purchase_cost_unit=np.exp(8.5133 + 0.9847*lnS-0.0561*(lnS**2))*(qs.CEPCI/567) 
            purchase_cost=number_dryers*(purchase_cost_unit)
        
        self.baseline_purchase_costs['Dryer'] = purchase_cost  
        self.F_D['Dryer']=self.F_P['Dryer']=self.F_M['Dryer']=1  
        self.F_BM['Grinder']=2.06

