
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.CfB_functions import formation_energy

cc = ComponentContribution.init()

listKeggMolecules = ['C00015']


print formation_energy(listKeggMolecules,cc, 8.25, 0.25)

