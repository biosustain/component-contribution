import math
import scipy.stats as st
from component_contribution.molecule import Molecule
from component_contribution.compound_cacher import KeggCompoundCacher
from component_contribution.CfB_functions import process_input_mets, convert2any
from pdbremix import pdbatoms
from pdbremix import asa
from component_contribution.thermodynamic_constants import default_RT
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.all_model import AllModel
from component_contribution.CfB_functions import add_thermo_comp_info
import os

# Physiological conditions
pH = 7.2
I = 0.1        # ionic strength (M)
T = 298.15     # temperature    (K)

# Confidence level
conf = 0.95

# INPUT Concentrations as file or string of MOL files
input_ids = ['h2o_c','xyl3_c','xyl__D_c','3oddecACP_h','h_h','nadph_h','3hddecACP_h','nadp_h','glyc3p_c','palmACP_c',\
             '1hdecg3p_c','ACP_c','h2o_c','C00018.sdf','C00009.sdf','C00250.sdf','atp_c','h2o_c','adp_c','pi_c','atp_c',\
             'pydx_c','adp_c','h_c','pydx5p_c']

input_ids = list(set(input_ids))
input_ids = [i for i in input_ids if not i.startswith('h_') and not i.startswith('h2o_')]

ccache = KeggCompoundCacher()

all_this = process_input_mets(input_ids,ccache)


# create one dictionary that contains all of the input_ids:compound instance
kegg_comps={}
for x in all_this['to_kegg_dict']:
    cid_id = all_this['to_kegg_dict'][x]
    kegg_comps[x] = ccache.get_compound(cid_id)
non_kegg_comps = all_this['non_kegg_refs']
all_comps = non_kegg_comps.copy()
all_comps.update(kegg_comps)

# calc charge and NPSA, and [] for all comps
countcharge = {}
NPSA = {}
concentration = {}

for input_id,comp in zip(all_comps.keys(), all_comps.values()):
    print input_id

    if comp == None:
        countcharge[input_id] = None
        NPSA[input_id] = None
        concentration[input_id] = None
        continue

    # calculate charge

    # charge calculations work better from smiles
    if comp.smiles_pH7:
        smiles_pH7 = comp.smiles_pH7
        mol = Molecule.FromSmiles(smiles_pH7)
    else:
        countcharge[input_id] = None
        NPSA[input_id] = None
        concentration[input_id] = None
        continue

    charge_list = mol.GetAtomCharges()
    charge_count = sum(x != 0 for x in charge_list)
    countcharge[input_id] = charge_count



    # calc NPSA

    # NPSA calculations work better from mol files
    molfile = comp.molfile
    mol = Molecule.FromMol(molfile)

    # create the pdb file. in the futre can cache this
    test = Molecule._ToFormat(mol.obmol,'pdb')
    file_name = "../examples/pdb_files/" + input_id +".pdb" # have to decide what we're going to save all the files as
    with open(file_name, "w") as text_file:
        text_file.write(test)

    # load pdb create " soup".     A Soup is essentially a collection of atoms, which we can grab by:
    soup = pdbatoms.Soup(file_name)
    atoms = soup.atoms()

    # asa - calculates the accessible surface - area of every atom in list of atom, with respect to the other atoms. which  assigns the asa to to each atom.asa
    pdbatoms.add_radii(atoms) # this calculates the radius of each atom.
    areas = asa.calculate_asa(atoms, 1.4) # 1.4 is the probe i.e. the size of water, changing this can change results ~ 25%
    total_area = sum(areas)

    # get atom neighbors
    adj = asa.adjacency_list(atoms, 1.8)
    adj = [[atoms[c].element for c in at] for at in adj]

    # get polar surface area, i.e. the area contributed by polar atoms only (oxygen, nitrogen and the hydrogen atoms attached to them
    polar_area=0
    for a, area, adj1 in zip(atoms, areas,adj):
        print a, area
        if a.element in ['O','N']:       # a.bfactor = area
            polar_area += area
        if a.element=='H' and  any([i in ['O','N'] for i in adj1]):
            polar_area += area

    NPSA1 = total_area - polar_area
    NPSA[input_id] = NPSA1

    conc = math.exp(charge_count*1.0425-NPSA1*0.0272)/1000

    concentration[input_id] = conc

print concentration, NPSA, countcharge





#d['InChI'] = comp.inchi
# try:
#     mol = Molecule.FromInChI(str(comp.inchi))
#     d['mass'] = mol.GetExactMass()
#     d['formula'] = mol.GetFormula()


# # Critical value for margin of error
# Z = st.norm.ppf(1-(1-conf)/2)
#
# # Print results
# for i, rxn_str in enumerate(reaction_strings):
#     print '-' * 60
#     print rxn_str.strip()
#     # Change in Gibbs free energy in standard conditions
#     print "dG0  = %8.1f +- %5.1f kj/mol" % (model.dG0[i, 0], dG0_std[i, 0] * Z)
#     # Change in Gibbs free energy at a particular pH and ionic strength
#     print "dG'0 = %8.1f +- %5.1f kj/mol" % (dG0_prime[i, 0], dG0_std[i, 0] * Z)
#     # Change in Gibbs free energy accounting for reactant concentration
#     print "dG'm = %8.1f +- %5.1f kj/mol" % (dGm_prime[i, 0], dG0_std[i, 0] * Z)