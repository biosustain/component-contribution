import os

import numpy as np
import scipy.stats as st

from component_contribution.core.metabolic_model import MetabolicModel
from component_contribution.prediction_model.model import ComponentContributionModel
from component_contribution.thermodynamics.constants import STANDARD_RT

# Physiological conditions
pH = 7.2
I = 0.1        # ionic strength (M)
T = 298.15     # temperature    (K)
conc_mM = 1.0  # comcentration  (mM)

# Confidence level
conf = 0.95

# INPUT Reactions as file or string of KEGG IDs
example_path = os.path.dirname(os.path.realpath(__file__))
REACTION_FNAME = os.path.join(example_path, 'wolf_reactions.txt')
with open(REACTION_FNAME, 'r') as fp:
    reaction_strings = fp.readlines()

# Parse the reaction strings into a "KEGG" prediction_model
model = MetabolicModel.from_kegg_formulas(reaction_strings, raise_exception=True)

# Train component-contribution predictive prediction_model
cc = ComponentContributionModel.init()

# Add the component-contribution estimation for reaction's untransformed dG0
# (i.e. using major MS at pH 7 for each of the reactants)
model.add_thermo(cc)

# Change in Gibbs free energy at a particular pH and ionic strength
dG0_prime, dG0_std, sqrt_Sigma = model.get_transformed_dG0(pH,I,T)

# Concentration matrix (M)
conc_M = 1e-3 * np.matrix(np.full((len(model.compound_ids), 1), conc_mM))
# Set Water "concentration" to 1 M
if 'C00001' in model.compound_ids:
    conc_M[model.compound_ids.index('C00001'), 0] = 1.0

# Change in Gibbs free energy accounting for reactant concentration
dGm_prime = dG0_prime + STANDARD_RT * model.s_matrix.T * np.log(conc_M)

# Critical value for margin of error
Z = st.norm.ppf(1-(1-conf)/2)

# Print results
for i, rxn_str in enumerate(reaction_strings):
    print('-' * 60)
    print(rxn_str.strip())
    # Change in Gibbs free energy in standard conditions
    print("dG0  = %8.1f +- %5.1f kj/mol" % (model.dG0[i, 0], dG0_std[i, 0] * Z))
    # Change in Gibbs free energy at a particular pH and ionic strength
    print("dG'0 = %8.1f +- %5.1f kj/mol" % (dG0_prime[i, 0], dG0_std[i, 0] * Z))
    # Change in Gibbs free energy accounting for reactant concentration
    print("dG'm = %8.1f +- %5.1f kj/mol" % (dGm_prime[i, 0], dG0_std[i, 0] * Z))