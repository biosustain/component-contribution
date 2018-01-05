import os

import numpy as np
import scipy.stats as st

from CfB_functions import add_thermo_comp_info
from component_contribution.core.metabolic_model import MetabolicModel
from component_contribution.prediction_model.model import ComponentContributionModel
from component_contribution.thermodynamics.constants import STANDARD_RT

# Physiological conditions
pH = 7.0
I = 0.1        # ionic strength (M)
T = 298.15     # temperature    (K)

# Confidence level
conf = 0.95

# INPUT Reactions as file or string of MOL files
example_path = os.path.dirname(os.path.realpath(__file__))

with open('../examples/example_bigg_id', 'r') as fp:
    reaction_strings = fp.readlines()

# Parse the reaction strings into a "KEGG" prediction_model
model = MetabolicModel.from_formulas(reaction_strings, 'bigg', raise_exception=True)

# Train component-contribution predictive prediction_model
cc = ComponentContributionModel.init()

# Add the component-contribution estimation for reaction's untransformed dG0
# (i.e. using major MS at pH 7 for each of the reactants)
# prediction_model.add_thermo(cc)
add_thermo_comp_info(model, cc)

# Change in Gibbs free energy at a particular pH and ionic strength
dG0_prime, dG0_std, sqrt_Sigma = model.get_transformed_dG0(pH, I, T)

# use 1 mM default Concentration
conc_M =  np.matrix(np.full((len(model.input_metabolites), 1), 1e-3))

# Set Water "concentration" to 1 M
if 'C00001' in model.compound_ids:
    for met, comp in model.comp_data['met_to_comp_dict'].iteritems():
        if comp == 'C00001':
            print(met)
            conc_M[ model.input_metabolites.index(met), 0] = 1.0


# Change in Gibbs free energy accounting for reactant concentration
dGm_prime = dG0_prime + STANDARD_RT * model.s_matrix_mets.T * np.log(conc_M)

# Critical value for margin of error
Z = st.norm.ppf(1-(1-conf)/2)

# Print results


file_name = "../examples/output_22.txt" # have to decide what we're going to save all the files as
with open(file_name, "w") as text_file:
    for i, rxn_str in enumerate(reaction_strings):
        print('-' * 60)
        print(rxn_str.strip())
        # Change in Gibbs free energy in standard conditions
        print("dG0  = %8.1f +- %5.1f kj/mol" % (model.dG0[i, 0], dG0_std[i, 0] * Z))
        # Change in Gibbs free energy at a particular pH and ionic strength
        print("dG'0 = %8.1f +- %5.1f kj/mol" % (dG0_prime[i, 0], dG0_std[i, 0] * Z))
        # Change in Gibbs free energy accounting for reactant concentration
        print("dG'm = %8.1f +- %5.1f kj/mol" % (dGm_prime[i, 0], dG0_std[i, 0] * Z))

        text_file.write(rxn_str.strip() + '\t' +str(model.dG0[i, 0])+ '\t'+ str(dG0_prime[i, 0]) + '\t' + str(dG0_std[i, 0]*Z)+'\n')