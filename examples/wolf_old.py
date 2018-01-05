import os

import numpy as np

from component_contribution.core.metabolic_model import MetabolicModel
from component_contribution.prediction_model.model import ComponentContributionModel
from component_contribution.thermodynamics.constants import STANDARD_RT

example_path = os.path.dirname(os.path.realpath(__file__))
REACTION_FNAME = os.path.join(example_path, 'wolf_reactions.txt')

cc = ComponentContributionModel.init()
with open(REACTION_FNAME, 'r') as fp:
    reaction_strings = fp.readlines()
model = MetabolicModel.from_kegg_formulas(reaction_strings, raise_exception=True)

model.add_thermo(cc)
dg0_prime, dg0_std, sqrt_sigma = model.get_transformed_dg0(7.0, 0.1, 298.15)

mM_conc = 1e-3 * np.matrix(np.ones((len(model.compound_ids), 1)))
if 'C00001' in model.compound_ids:
    mM_conc[model.compound_ids.index('C00001'), 0] = 1.0
dGm_prime = dg0_prime + STANDARD_RT * model.s_matrix.T * np.log(mM_conc)

for i, r in enumerate(reaction_strings):
    print('-'*50)
    print(r.strip())
    print("dG0  = %8.1f +- %5.1f" % (model.dg0[i, 0], dg0_std[i, 0] * 1.96))
    print("dG'0 = %8.1f +- %5.1f" % (dg0_prime[i, 0], dg0_std[i, 0] * 1.96))
    print("dG'm = %8.1f +- %5.1f" % (dGm_prime[i, 0], dg0_std[i, 0] * 1.96))
