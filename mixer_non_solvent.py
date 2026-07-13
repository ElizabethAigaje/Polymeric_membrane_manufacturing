
"""
LC systems analysis to identify green polymer materials and methods for membrane manufacturing

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 12/5/2025
"""
from thermosteam._graphics import mixer_graphics
import flexsolve as flx
import biosteam as bst
import numpy as np
from typing import Optional
import qsdsan as qs

__all__ = ('Mixer', 'SteamMixer', 'FakeMixer', 'MockMixer')



class Mixer_non_solvent(qs.SanUnit):
    """ 
    Dummy mixer to calculate non-solvent mass flow as function of the plant capacity. 
    
    Parameters
    ----------
    ins : 
        Inlet non-solvent, (the default is water).
    outs : 
        Outlet non-solvent based on capacity.
    nonsolvent_per_membrane :
        Mass of non-solvent per membrane area, [kg/m^2] (the default is 150 kg/m^2 [1])
    params:
        Variable that stores system-level parameters, [-] (see class 'MembraneParams' in the system-membrane-manufacutring.py file).
    rigorous :
        Whether to perform vapor-liquid equilibrium.
    
    Notes
    -----
    When streams at different pressures are mixed, BioSTEAM assumes valves 
    reduce the pressure of the streams being mixed to prevent backflow 
    (pressure needs to decrease in the direction of flow according to 
    Bernoulli's principle). The outlet pressure will be the minimum pressure
    of all inlet streams.
    The amount of water required for coagulation was calculated based on reference [1], that indicates that wastewater generation in a 
    membrane manufacturing plant can be estimated to be 100-500 L of wastewater. Therefore, the midpoint is 300 L (300 kg) produced during 
    coagulation and rising. We assume half of the water produced here and the other half in the rinsing tank.

    References
    ----------
    [1] Razali, Ml, Kim, J., Attfield, M., Budd, P., Drioli, E., Lee, Y. and Szekely, G. Sustainable wastewater treatment and recycling in membrane manufacturing.
        Green Chem., 2015, 17, 5196.
    """
    
    _graphics = mixer_graphics
    _N_outs = 1
    _N_ins = 1
    _ins_size_is_fixed = False
    
    def _assert_compatible_property_package(self): 
        pass # Not necessary for mixing streams
    
    def _init(self, rigorous: Optional[bool]=False,
              conserve_phases: Optional[bool]=False, 
               nonsolvent_per_membranearea=150 , params=None):  
        self.rigorous = rigorous
        self.nonsolvent_per_membranearea=nonsolvent_per_membranearea
        self.conserve_phases = conserve_phases
        self.params=params
        
    def _run(self):
        p=self.params
        polymer_needed=(p.v_polymersol * p.d_polymersol *
                         p.polymer_fraction) / 1000.0 / 24.0  #kg/h of polymer
        non_solvent_needed = (polymer_needed/p.polymermass_to_membranearea)*self.nonsolvent_per_membranearea     #kg/h of non-solvent

    
    # obtain inlet stream object safely (handles N_ins == 1 or >1)
        inlet = self.ins[0] 

        # clear the inlet
        inlet.empty()

        # Set the mass flow correctly
        # You MUST give a dict with component ID
        inlet.imass['H2O'] = non_solvent_needed

        # get outlet and clear
        outlet = self.outs[0]
        outlet.empty()

    # mix
        outlet.mix_from([inlet], vle=self.rigorous,
                conserve_phases=getattr(self, 'conserve_phases', None))

        V = outlet.vapor_fraction
        if V == 0:
            self._B = 0
        elif V == 1:
            self._B = np.inf
        else:
            self._B = V / (1 - V)
    
    
    def _get_energy_departure_coefficient(self, stream):
        if stream.phases == ('g', 'l'):
            vapor, liquid = stream
            if vapor.isempty():
                with liquid.temporary_phase('g'): coeff = liquid.H
            else:
                coeff = -vapor.h * liquid.F_mol
        else:
            coeff = -stream.C
        return (self, coeff)
    
    def _create_energy_departure_equations(self):
        # Ll: C1dT1 - Ce2*dT2 - Cr0*dT0 - hv2*L2*dB2 = Q1 - H_out + H_in
        # gl: hV1*L1*dB1 - hv2*L2*dB2 - Ce2*dT2 - Cr0*dT0 = Q1 + H_in - H_out
        outlet = self.outs[0]
        phases = outlet.phases
        if phases == ('g', 'l'):
            vapor, liquid = outlet
            coeff = {}
            if vapor.isempty():
                with liquid.temporary_phase('g'): coeff[self] = liquid.H
            else:
                coeff[self] = vapor.h * liquid.F_mol
        else:
            coeff = {self: outlet.C}
        for i in self.ins: i._update_energy_departure_coefficient(coeff)
        return [(coeff, self.H_in - self.H_out)]
    
    def _create_material_balance_equations(self, composition_sensitive):
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        outlet, = self.outs
        if len(outlet) == 1:
            ones = np.ones(self.chemicals.size)
            minus_ones = -ones
            zeros = np.zeros(self.chemicals.size)
            
            # Overall flows
            eq_overall = {outlet: ones}
            for i in process_inlets: eq_overall[i] = minus_ones
            equations.append(
                (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
            )
        else:
            top, bottom = outlet
            ones = np.ones(self.chemicals.size)
            minus_ones = -ones
            zeros = np.zeros(self.chemicals.size)
            
            # Overall flows
            eq_overall = {}
            for i in outlet: 
                eq_overall[i] = ones
            for i in process_inlets:
                eq_overall[i] = minus_ones
            equations.append(
                (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
            )
            
            # Top to bottom flows
            B = self._B
            eq_outs = {}
            if B == np.inf:
                eq_outs[bottom] = ones
            elif B == 0:
                eq_outs[top] = ones
            else:
                bp = outlet.bubble_point_at_P()
                outlet.T = bp.T
                S = bp.K * B
                eq_outs[top] = ones
                eq_outs[bottom] = -S
            equations.append(
                (eq_outs, zeros)
            )
        return equations
    
    def _update_energy_variable(self, departure):
        phases = self.outs[0].phases
        if phases == ('g', 'l'):
            self._B += departure
        else:
            self.outs[0].T += departure

    def _update_nonlinearities(self): pass