"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 03/24/2026
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom

class Dilution_tank(qs.sanunits.MixTank): 
    """ Creates an agitated tank to dilute wastewater before its discharge to off-site treatment,

    Parameters
    ----------
    ID : str
        ID of the dissolution tank.
    ins : Iterable(stream)
        sodiumhydroxide, water.
    outs: Iterable(stream)
        dissolution.
    COD_target: float
        Max COD permisible by municipal wastewater facilities, [mg/L] (the default is 425.36 [mg/L] [1]).
    tau : float
        Dissolution time, [hr]
    kW_m3 : float
        Power required for mixing (the default for blending, mixing is 0.07 [2]).

    Notes
    -----
    The EPA priority chemical list does not specify a limit concentration discharge for NMP. Morevoer, industry practices report
    dilution as a common on-site 'wastewater treatment' practice  in which the parameters taken in account are PH (6-9), COD, BOD, etc.
    Our wastewater has a neutral PH, so neutralizaiton does not taken place. We do attempt to reduce the COD to the average value reported in ref [1],
    before its discharge to off-site treatment.
    In ref [3], they indicate there is a limit of 0.1% for solvent residues stablished by the Restriction of Chemical regulations of the 
    European Union countries. For the baseline, the final % we get for the non-reclycing scenario is 0.02% with the water calculated to reduce 
    the COD, so we are even below that concentration.
    COD may be a more relevant parameter as NMP is not in the list or controled in wastewater treatment facilities in the USA.

    References
    ---------
    [1] Ecoinvent database. Activity: Treatment of wastwater, average wastewater treatment [ROW]
    [2] Sinnott, R. K., & Towler, G. P. (2020). Chemical engineering design (6th edition). Elsevier, pp. 625.
    [3] Li.N, et al., Biodegradation of N-methyl-2-pyrrolidone (NMP) in wastewater: A review of current knowledge and future perspectives.
        Journal of Cleaner Production 486 (2025) 144452.
        
    """
    
    _N_ins=2
    _N_outs=1


    def __init__(self, ID, ins=(), outs=(), Thermo=None, 
                 COD_target=425.34, tau=(), kW_per_m3=2):   #max solubility at 120 C for LPDE in solvent, power required to mix slurry suspensions 
        super().__init__(ID, ins, outs, Thermo) 
        
        self.COD_target=COD_target
        self.tau=tau
        self.kW_per_m3=kW_per_m3

    def _run(self):
        """Calculate material flow distribution and phase conditions."""
        wastewater, water = self.ins           
        diluted_wastewater, = self.outs               
        
        # Calculate the initial COD 
       
        wastewater_nmp= wastewater.imass['NMP']
        wastewater_water=wastewater.imass['H2O']
        wastewater_etoh=wastewater.imass['EtOH']  #kg/h

        nmp_fraction= wastewater_nmp/(wastewater_nmp+wastewater_water)
        etoh_concentration = wastewater_etoh /wastewater.F_vol/(1000)  #mg/L = kg/h/(m3/h)

        NMP_COD = 1908776 #mg/L from simulation of NMP stream in jupiternotebook using qsdsan
        
        #COD from ethanol
    
        etOH_CODload =2.09*wastewater_etoh # kg O2, 2.09 by chat GPT 
        etOH_COD=(etOH_CODload/wastewater.F_vol)*1000  #kg/m^3 --> mg/L
       
        wastewater_COD = NMP_COD*nmp_fraction + etOH_COD

        # Calculate the water needed to get to the target COD
        if wastewater_COD > self.COD_target:   #The if is needed specially for the recovery cases
            water_needed=wastewater.F_mass*((wastewater_COD/self.COD_target)-1)   #kg/h
            water.imass['H2O']=water_needed

        else:
            water.imass['H2O'] == 0
            
        
        diluted_wastewater.copy_like(wastewater)
            #self.outs[0].copy_like(self.ins[0])
        
        diluted_wastewater.imass['H2O']= wastewater.imass['H2O']+water.imass['H2O']

        # print(self.ID, "RUNNING")
        # print("h2o kg/h =", water_needed)
        # print("etoh COD =", etOH_COD)
        # print('nmp_cod=', NMP_COD*nmp_fraction)
        

    def _design(self):
        """Inherits design calculation from MixTank."""
        qs.sanunits.MixTank._design(self)
        #Accounts for the number of tanks needed to work in paralell to cover continous process

    def _cost(self):
        """Inherits cost calculation from MixTank."""
        qs.sanunits.MixTank._cost(self)
       
         
         



   
