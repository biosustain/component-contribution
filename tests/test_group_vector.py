from component_contribution.core import group_vector

atp_inchi = 'InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29' \
            '(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1/f/h18-' \
            '19,21,23H,11H2'
atp_smiles_ph7 = 'NC1=C2N=CN([C@@H]3O[C@H](COP([O-])(=O)OP([O-])(=O)OP(O)([O-])=O)[C@@H](O)[C@H]3O)C2=NC=N1'
group_dict = {43: 1, 50: 1, 72: 1, 74: 3, 103: 1, 119: 2, 127: 1, 147: 2, 148: 2, 149: 1, 157: 2, 160: 1, 162: 1}


def test_decomposer():
    groups_data = group_vector.init_groups_data()
    group_decomposer = group_vector.GroupDecomposer(groups_data)
    groupvec1 = group_vector.inchi_to_group_vector(group_decomposer, atp_inchi)
    groupvec2 = group_vector.smiles_to_group_vector(group_decomposer, atp_smiles_ph7)
    print(groupvec1)
    print(groupvec2)
    for group_ind, group_count in enumerate(groupvec1.flatten()):
        assert(group_dict.get(group_ind, 0) == group_count)
