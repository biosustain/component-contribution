# This script compares CfB's version of component contribution with equilibrator.
#
# up until line 44 this script only takes '../../validation/internal/kegg_reactions', a text file with all kegg reactions
# and transforms it to '../../validation/internal/kegg_to_met_reactions' which replaces metabolite names with bigg_ids
# the last line

# the last line actually calls the function to calculate reaction gibbs energy as in point_23_4


import json
import re

from CfB_functions import load_bigg_dict, output_reaction_energies

# # first load bigg database to translate. bigg dict is the downloaded bigg metabolite database
bigg_dict = load_bigg_dict(filename='../data/bigg_models_metabolites.tsv')
keggcomp_to_met_dict = {}

# map all ids I can to kegg (or whatever internal ID with Compound info).
for met in bigg_dict.keys():
    try:
        # classifying met. first see if met matches the bigg database
        all_references_json = bigg_dict[met]
        all_references_readable = json.loads(all_references_json)

        if 'KEGG Compound' in all_references_readable: # we matched it to kegg met and will use component contribution's database for this guy
            kegg_reference = all_references_readable['KEGG Compound']
            kegg_id = kegg_reference[0]['id']
            keggcomp_to_met_dict[kegg_id] = met
    except:
        pass

    # get all the kegg reactions
with open('../../validation/internal/kegg_reactions', 'r') as fp:
    allreac_names = fp.readlines()
allreac_names = [x.strip() for x in allreac_names]

with open('../../validation/internal/kegg_to_met_reactions', 'w') as fp:
    for i, reac in enumerate(allreac_names):
        if i%200==0: print(i)
        temp_reac=reac
        temp_reac.replace('<==>', '+')
        pieces = re.split(' ', temp_reac)

        for piece in pieces:
            if piece in keggcomp_to_met_dict.keys():
                temp_reac = temp_reac.replace(piece, keggcomp_to_met_dict[piece])
        fp.write(temp_reac + '\n')

output_reaction_energies('../../validation/internal/kegg_to_met_reactions',7,  0.1,'../../validation/internal/kegg_to_met_reactions_output')
