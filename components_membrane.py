"""
File of the list of components

LC systems analysis to identify green polymer materials and methods
for membrane manufacturing.

Developed by: Elizabeth Aigaje-Espinosa
Chemical Engineering Department
The Pennsylvania State University
Advisor: Rui Shi

Last modified: 04/27/2026
"""


import qsdsan as qs

__all__ = ('create_components_membrane', )


def create_components_membrane():
    """
    Create the components needed for membrane manufacturing.

    Returns
    -------
    comps1 : qs.Components
        Components ready to be used in the simulation.
    """
    
#Polymer and bore solution
    H2O_chem=qs.Chemical('H2O')
    H2O = qs.Component.from_chemical(ID='H2O', chemical=H2O_chem, particle_size='Soluble', degradability='Undegradable', organic=False)

    N2_chem=qs.Chemical('N2')
    N2 = qs.Component.from_chemical(ID='N2', chemical=N2_chem, particle_size='Soluble', degradability='Undegradable', organic=False)

    NMP_chem=qs.Chemical('NMP')
    NMP = qs.Component.from_chemical(ID='NMP', chemical=NMP_chem, particle_size='Soluble', degradability='Undegradable', organic=True)
 

    
    polysulfone = qs.Component(ID='polysulfone',  phase='s', particle_size='Particulate', degradability='Undegradable', organic=True,  measured_as = 'COD', MW=60000)
    heptacontane_chem = qs.Chemical('heptacontane')

    #Surrogate for thermo models
    polysulfone.Tm=heptacontane_chem.Tm
    polysulfone.Tb=heptacontane_chem.Tb
    polysulfone.Pc=heptacontane_chem.Pc
    polysulfone.Tc=heptacontane_chem.Tc
    polysulfone.Vc=heptacontane_chem.Vc
    polysulfone.Hfus=heptacontane_chem.Hfus
    polysulfone.Hf=heptacontane_chem.Hf
    polysulfone.copy_models_from(heptacontane_chem)
    polysulfone.V.add_method(0.04838)   #60000/1240/1000  gr/mol/kg/m^3
    # polysulfone.Cn.add_model(82200)   #molar heat capacity 1.37 J/gk * 60000 g/mol
    polysulfone.Psat.add_method(0)
    #Consider essentially solid resins that do not fit assumption of the EOS/VLE thermodynamic models in QSDsan  

#Additives

    PVP_chem=qs.Chemical('PVP')
    PVP = qs.Component.from_chemical(ID='PVP', chemical=PVP_chem, particle_size='Soluble', degradability='Undegradable', organic=False)

    PEG_chem=qs.Chemical('PEG')
    PEG = qs.Component.from_chemical(ID='PEG', chemical=PEG_chem, particle_size='Soluble', degradability='Undegradable', organic=True)


    glycerol_chem=qs.Chemical('glycerol')
    glycerol = qs.Component.from_chemical(ID='glycerol', chemical=glycerol_chem, particle_size='Soluble', degradability='Undegradable', organic=True)

#Adhesives
#Epoxy resin Bisphenol A epoxy resin:  https://commonchemistry.cas.org/detail?ref=25068-38-6
    epoxy_chem=qs.Chemical('araldite 527')
    epoxy = qs.Component.from_chemical(ID='epoxy', chemical=epoxy_chem, particle_size='Soluble', degradability='Undegradable', organic=False)
    epoxy.Hfus=heptacontane_chem.Hfus
    epoxy.Hf=heptacontane_chem.Hf    #it really does not impact as it does not goes though any calculation
                                     # but it needs to be defined because if not it causes an error in the heat exchanger calculations

#Cellulose regeneration
    NaOH_chem=qs.Chemical('NaOH')
    NaOH = qs.Component.from_chemical(ID='NaOH', chemical=NaOH_chem, particle_size='Particulate', degradability='Undegradable', organic=False)

    EtOH_chem=qs.Chemical('ethanol')
    EtOH = qs.Component.from_chemical(ID='EtOH', chemical=EtOH_chem, particle_size='Soluble', degradability='Undegradable', organic=False)

    cmps1 = qs.Components([H2O, N2, NMP, polysulfone, PVP, PEG, glycerol, epoxy, NaOH, EtOH]) # *cmps_default
    qs.set_thermo(cmps1, cache=True)
   

    return cmps1
    
cmps1=create_components_membrane()
print(cmps1)
