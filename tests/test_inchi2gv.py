import sys, logging
sys.path.append('../python')
from compound import Compound
import inchi2gv
from compound_cacher import CompoundCacher
from molecule import Molecule

#logger = logging.getLogger('')
#logger.setLevel(logging.DEBUG)
ccache = CompoundCacher('../cache/compounds.json')
groups_data = inchi2gv.init_groups_data()
group_list = groups_data.GetGroupNames()
decomposer = inchi2gv.InChIDecomposer(groups_data)

patterns = ['c~[O;+0]', 'c~[O;+1]', 'c~[n;+1]~c', 'c~[n;+0]~c', 'c~[n;-1]~c']

for cid in [255, 1007]:
    comp = ccache.get_kegg_compound(cid)
    print "-"*50, '\nC%05d' % cid
    inchi = comp.inchi
    mol = Molecule.FromInChI(inchi)
    print mol.ToSmiles()
    
    print mol.FindSmarts("c~[n;+1]~c")
    
    try:
        groupvec = decomposer.inchi_to_groupvec(inchi)
        sys.stdout.write(str(groupvec) + '\n')
    except inchi2gv.GroupDecompositionError as e:
        sys.stderr.write(str(e) + '\n')
        sys.stderr.write(e.GetDebugTable())
    