from component_contribution.CfB_functions import output_reaction_energies

# this script calls the function "output_reaction_energies" that puts all of the previously developed functionality.
# inputs:
#  text file with all the reactions one would like to calculate gibbs reaction energy for
#  pH
#  Ionic strength
#  file name to output the resulting information. the written text file will include the columns:
#       reaction formula	model.dG0	dG0_prime	dGm_prime	dG0_std

# here I am calculating the the reaction energies for E.coli iAF1260 at different conditions. I then compare these results
# to the Hatzimanikatis thermo data in the same model in the script "point_25"

output_reaction_energies('../../validation/allreacs.txt', 7,   0.1, '../examples/pH7_I01_edited.txt')
output_reaction_energies('../../validation/allreacs.txt', 7,   0,   '../examples/pH7_I0_edited.txt')
output_reaction_energies('../../validation/allreacs.txt', 7.2, 0.1, '../examples/pH72_I01_edited.txt')
output_reaction_energies('../../validation/allreacs.txt', 7.2, 0,   '../examples/pH72_I0_edited.txt')