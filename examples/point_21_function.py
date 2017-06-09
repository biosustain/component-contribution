# This script is an example of how we would calculate the gibbs formation of a compound. the input file contains
# reactions with single metabolites on the product side.

# this method is not recommended by Elad Noor and Ron Milo.


import numpy as np
import scipy.stats as st
from component_contribution.thermodynamic_constants import default_RT
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.all_model import AllModel
from component_contribution.CfB_functions import add_thermo_comp_info, calc_concentration
import os

# Physiological conditions
pH = 7.2
I = 0.1        # ionic strength (M)
T = 298.15     # temperature    (K)
conc_mM = 1.0  # comcentration  (mM)

# Confidence level
conf = 0.95

# INPUT Reactions as file or string of MOL files
example_path = os.path.dirname(os.path.realpath(__file__))
with open('../examples/example_bigg_mets', 'r') as fp:

    reaction_strings = fp.readlines()

# Parse the reaction strings into a "KEGG" model
model = AllModel.from_formulas(reaction_strings, 'bigg', raise_exception=True)

# Train component-contribution predictive model
cc = ComponentContribution.init()

# Add the component-contribution estimation for reaction's untransformed dG0
# (i.e. using major MS at pH 7 for each of the reactants)
# model.add_thermo(cc)
add_thermo_comp_info(model, cc)

# Change in Gibbs free energy at a particular pH and ionic strength
dG0_prime, dG0_std, sqrt_Sigma = model.get_transformed_dG0(pH,I,T)

# Critical value for margin of error
Z = st.norm.ppf(1-(1-conf)/2)

# Print results
for i, rxn_str in enumerate(reaction_strings):
    print '-' * 60
    print rxn_str.strip()
    #  Gibbs formation energy in standard conditions
    print "dG0  = %8.1f +- %5.1f kj/mol" % (model.dG0[i, 0], dG0_std[i, 0] * Z)
    # Gibbs formation energy at a particular pH and ionic strength
    print "dG'0 = %8.1f +- %5.1f kj/mol" % (dG0_prime[i, 0], dG0_std[i, 0] * Z)
