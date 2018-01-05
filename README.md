README for CfB Component contribution.

This package is a fork of https://github.com/eladnoor/component-contribution

It is modified to take as inputs biochemical reactions written in bigg_ids and molfiles. It includes example scripts covering the points listed below:

point_21_function: this script demostrates accomplishment of Gibbs formation calculation

point_22_function: this script demostrates accomplishment of Gibbs reaction calculation
	This script exemplifies how the reaction gibbs energy can be calculated for any biochemical reaction that is made up of bigg_ids and mol files. when a metabolite is a mol file, the metabolite id in the reaction is its file name and the mol file must be in examples/molfiles (this can be edited in the function process_input_mets)


point_23_4_function: this script estimates the concentration of the metabolites to improve reaction energy predictions.
	This script calls the function "output_reaction_energies" which takes an input text file with all the reactions one would like to calculate gibbs reaction energy for, the pH, Ionic strength, and the file name to output the resulting information. the written text file will include the columns: reaction formula	model.dG0	dG0_prime	dGm_prime	dG0_std

	Inside the "output_reaction_energies" is the calc_concentration function. If the first input to this function is True, it will calculate the concetration based on the metabolite structural properties like in Bar-Even 2011. If False, it will output 1 millimolar for almost all metabolites (see code for exceptions)

	This script runs output_reaction_energies four times on the whole iAF1260 model in four different conditions to then be compared with the  Thermo data in that model. output files are generated in examples/

point_25: this script shows a validation of the method
	This script runs the validation against Feist et. al. data. It takes precalculated data from "point_23_4_function" and just does some filtering to properly match with Feist data. RESULTS ARE THEN MANUALLY COMPILED IN THE SPREADSHEET "compile_point25" in the examples folder

point_25internal: this script shows an internal validation of the method
	This script compares CfB's version of component contribution with equilibrator. Up until line 44 this script only takes '../../validation/internal/kegg_reactions', a text file with all kegg reactions and transforms it to '../../validation/internal/kegg_to_met_reactions' which replaces metabolite names with bigg_ids.
	The last line actually calls the function to calculate reaction gibbs energy as in point_23_4. The resulting data was compiled in the spreadsheet "equilibrator_reaction_energies" in validation/internal/


==========================

### Requirements:
* python == 2.7
* numpy >= 1.6.2
* scipy >= 0.14
* oct2py >= 3.1 (dependent on Octave)
* Open-Babel >= 2.3.1
* ChemAxon's Marvin >= 5.11

in addition need the following
Pandas: to run script point_25
chebi_web: to translate bigg to kegg via chebi. this is commented out for now.
