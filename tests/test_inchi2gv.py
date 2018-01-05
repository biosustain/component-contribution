import sys

from component_contribution.core.group_vector import init_groups_data, inchi_to_group_vector
from component_contribution.core.group_vector import GroupDecomposer, GroupDecompositionError
from component_contribution.compound_cache import CompoundCache
from component_contribution.core.molecule import Molecule


ccache = CompoundCache('../cache/compounds.json')

groups_data = init_groups_data()
group_list = groups_data.group_names()
group_names = groups_data.group_names()
decomposer = GroupDecomposer(groups_data)

# test the decomposition of ATP into groups
atp_inchi = ccache.get_compound('C00002').inchi
group_def = inchi_to_group_vector(decomposer, atp_inchi)
for j, group_name in enumerate(group_names):
    if group_def[j] != 0:
        print(group_name, ' x %d' % group_def[j])


patterns = ['c~[O;+0]', 'c~[O;+1]', 'c~[n;+1]~c', 'c~[n;+0]~c', 'c~[n;-1]~c']

for cid in ['C00255', 'C01007']:
    comp = ccache.get_compound(cid)
    print("-"*50, '\n%s' % cid)
    inchi = comp.inchi
    mol = Molecule.from_inchi(inchi)
    print(mol.ToSmiles())
    
    print(mol.FindSmarts("c~[n;+1]~c"))
    
    try:
        groupvec = inchi_to_group_vector(decomposer, inchi)
        sys.stdout.write(str(groupvec) + '\n')
    except GroupDecompositionError as e:
        sys.stderr.write(str(e) + '\n')
        sys.stderr.write(e.debug_table())
    