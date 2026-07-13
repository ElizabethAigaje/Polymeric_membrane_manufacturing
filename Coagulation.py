"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 12/08/2025
"""

import qsdsan as qs
import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year
from biosteam.units.decorators import cost
from qsdsan.utils import auom
import math


class Coagulation_bath(qs.SanUnit):
    """ Creates a coagulation bath system for fiber precipitation (NIPS).

    Parameters
    ----------
    ID : str
        ID of extruder.
    ins : Iterable(stream)
        Inlet.
    outs: Iterable(stream)
        Outlet.
    solvent_out :
        Mass fraction of solvent difussing out due to induced phase separation, (the default is 0.7)


    Notes  
    -----
    Costs and energy demand for this stage, including winding, are accounted in the extruder cost.

    References
    ----------
  
    """
      

    def __init__(self, ID, ins=(), outs=(), thermo=None, solvent_out=0.7):
        super().__init__(ID, ins, outs,thermo)  
    
        self.solvent_out=solvent_out
       
    _N_ins=2
    _N_outs=2

    def _run(self):

        lumen_solution, nonsolvent_solution = self.ins           
        lumen_out, nonsolvent_and_solvent= self.outs               

        lumen_out.imass['polysulfone']=lumen_solution.imass['polysulfone']
        lumen_out.imass['PEG']=lumen_solution.imass['PEG']
        lumen_out.imass['H2O']= lumen_solution.imass['H2O']
        lumen_out.imass['NMP']= (1-self.solvent_out)*(lumen_solution.imass['NMP']+nonsolvent_solution.imass['NMP'])

        nonsolvent_and_solvent.imass['H2O']=nonsolvent_solution.imass['H2O']
        nonsolvent_and_solvent.imass['NMP']=self.solvent_out*(lumen_solution.imass['NMP']+nonsolvent_solution.imass['NMP'])
 
        lumen_out.phase ='l'
        nonsolvent_and_solvent.phase ='l'

    def _design(self):

        # The number of coagulation baths is included in the extrusion system
        pass
    def _cost(self):
        
        #The cost of coagulation baths is included in the extrusion system
        pass
