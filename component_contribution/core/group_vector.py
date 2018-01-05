import os
import sys

import csv
import json
import itertools

import logging

import numpy as np
from component_contribution.core.molecule import Molecule
from component_contribution.exceptions import OpenBabelError
from component_contribution.exceptions import GroupDecompositionError, GroupsDataError, MalformedGroupDefinitionError

logger = logging.getLogger(__name__)

GROUP_CSV = os.path.join(os.path.dirname(__file__), "groups.csv")


class GroupVector(list):
    """A vector of groups."""

    def __init__(self, groups_data, iterable=None):
        """Construct a vector.
        
        Arguments
        ---------
        groups_data : ???
            Data about all the groups.
        iterable : iterable
            Data to load into the vector.
        """
        super(GroupVector, self).__init__([])

        self.groups_data = groups_data
        
        if iterable is not None:
            self.extend(iterable)
        else:
            for _ in range(len(self.groups_data.all_group_names)):
                self.append(0)
    
    def __str__(self):
        """Return a sparse string representation of this group vector."""
        group_strs = []
        gv_flat = self.flatten()
        for i, name in enumerate(self.groups_data.group_names):
            if gv_flat[i]:
                group_strs.append('%s x %d' % (name, gv_flat[i]))
        return " | ".join(group_strs)
    
    def __iadd__(self, other):
        for i in range(len(self.groups_data.all_group_names)):
            self[i] += other[i]
        return self

    def __isub__(self, other):
        for i in range(len(self.groups_data.all_group_names)):
            self[i] -= other[i]
        return self
            
    def __add__(self, other):
        result = GroupVector(self.groups_data)
        for i in range(len(self.groups_data.all_group_names)):
            result[i] = self[i] + other[i]
        return result

    def __sub__(self, other):
        result = GroupVector(self.groups_data)
        for i in range(len(self.groups_data.all_group_names)):
            result[i] = self[i] - other[i]
        return result
    
    def __eq__(self, other):
        for i in range(len(self.groups_data.all_group_names)):
            if self[i] != other[i]:
                return False
        return True
    
    def __nonzero__(self):
        for i in range(len(self.groups_data.all_group_names)):
            if self[i] != 0:
                return True
        return False
    
    def __mul__(self, other):
        try:
            c = float(other)
            return GroupVector(self.groups_data, [x*c for x in self])
        except ValueError:
            raise ValueError("A GroupVector can only be multiplied by a scalar"
                             ", given " + str(other))

    @property
    def net_charge(self):
        """Returns the net charge."""
        return int(np.dot(self, self.groups_data.all_group_charges))

    @property
    def number_of_protons(self):
        """Returns the number of protons."""
        return int(np.dot(self, self.groups_data.all_group_hydrogens))

    @property
    def number_of_magnesium_ions(self):
        """Returns the number of Mg2+ ions."""
        return int(np.dot(self, self.groups_data.all_group_mgs))

    def remove_epsilon_values(self, epsilon=1e-10):
        for i in range(len(self)):
            if abs(self[i]) < epsilon:
                self[i] = 0

    def to_json(self):
        return json.dumps(dict([(i, x) for (i, x) in enumerate(self) if x != 0]))
    
    @staticmethod
    def from_json(groups_data, s):
        v = [0] * groups_data.count()
        for i, x in json.loads(s).iteritems():
            v[int(i)] = x
        return GroupVector(groups_data, v)
    
    def flatten(self):
        if not self.groups_data.transformed:
            return tuple(self)
        
        # map all pseudoisomeric group indices to Biochemical group indices (which are fewer)
        # use the names of each group and ignore the nH, z and nMg.
        biochemical_group_names = self.groups_data.group_names()
        biochemical_vector = [0] * len(biochemical_group_names)
        for i, x in enumerate(self):
            group_name = self.groups_data.all_groups[i].name
            new_index = biochemical_group_names.index(group_name)
            biochemical_vector[new_index] += x
        return tuple(biochemical_vector)        

    def to_array(self):
        return np.matrix(self.flatten())


class _AllAtomsSet(object):
    """A set containing all the atoms: used for focal atoms sets."""
    
    def __contains__(self, elt):
        return True


class FocalSet(object):
    def __init__(self, focal_atoms_str):
        if not focal_atoms_str:
            raise ValueError(
                'You must supply a non-empty focal atom string.'
                ' You may use "None" or "All" in the obvious fashion.')
        
        self.str = focal_atoms_str
        self.focal_atoms_set = None
        prepped_str = self.str.strip().lower()
        
        if prepped_str == 'all':
            self.focal_atoms_set = _AllAtomsSet()
        elif prepped_str == 'none':
            self.focal_atoms_set = set()
        else:
            self.focal_atoms_set = set([int(c) for c in self.str.split('|')])
    
    def __str__(self):
        return self.str
    
    def __contains__(self, elt):
        return self.focal_atoms_set.__contains__(elt)


class Group(object):
    """Representation of a single group."""
    
    def __init__(self, group_id, name, number_of_hydrogens, charge, number_of_magnesium_ions,
                 smarts=None, focal_atoms=None):
        self.id = group_id
        self.name = name
        self.number_of_protons = number_of_hydrogens
        self.charge = charge
        self.number_of_magnesium_ions = number_of_magnesium_ions
        self.smarts = smarts
        self.focal_atoms = focal_atoms

    def _is_hydrocarbon_group(self):
        return self.name.startswith('*Hc')

    def _is_sugar_group(self):
        return self.name.startswith('*Su')
    
    def _is_aromatic_ring_group(self):
        return self.name.startswith('*Ar')
    
    def _is_heteroaromatic_ring_group(self):
        return self.name.startswith('*Har')

    def is_phosphate(self):
        return self.name.startswith('*P')
    
    def ignore_charges(self):
        # (I)gnore charges
        return self.name[2] == 'I'
    
    def charge_sensitive(self):
        # (C)harge sensitive
        return self.name[2] == 'C'
    
    def is_coded_correction(self):
        """Returns True if this is a correction for which hand-written code.
           must be executed.
        """
        return (self._is_hydrocarbon_group() or
                self._is_aromatic_ring_group() or
                self._is_heteroaromatic_ring_group())

    @staticmethod
    def _is_hydrocarbon(mol):
        """Tests if a molecule is a simple hydrocarbon."""
        if mol.FindSmarts('[!C;!c]'):
            # If we find anything other than a carbon (w/ hydrogens)
            # then it's not a hydrocarbon.
            return 0
        return 1    

    @staticmethod
    def _count_aromatic_rings(mol):
        expressions = ['c1cccc1', 'c1ccccc1']
        count = 0
        for smarts_str in expressions:
            count += len(mol.FindSmarts(smarts_str))
        return count
    
    @staticmethod
    def _count_heteroaromatic_rings(mol):
        expressions = ['a1aaaa1', 'a1aaaaa1']
        count = 0
        all_atoms = mol.GetAtoms()
        for smarts_str in expressions:
            for match in mol.FindSmarts(smarts_str):
                atoms = set([all_atoms[i].atomicnum for i in match])
                atoms.discard(6)  # Ditch carbons
                if atoms:
                    count += 1
        return count

    def correction(self, mol):
        """Get the value of the correction for this molecule."""
        if self._is_hydrocarbon_group():
            return self._is_hydrocarbon(mol)
        elif self._is_aromatic_ring_group():
            return self._count_aromatic_rings(mol)
        elif self._is_heteroaromatic_ring_group():
            return self._count_heteroaromatic_rings(mol)
        
        raise TypeError('This group is not a correction.')
    
    def focal_set(self, nodes):
        """Get the set of focal atoms from the match.
        
        Args:
            nodes: the nodes matching this group.
        
        Returns:
            A set of focal atoms.
        """        
        focal_set = set()
        for i, node in enumerate(nodes):
            if i in self.focal_atoms:
                focal_set.add(node)            
        return focal_set
    
    def __str__(self):
        if self.number_of_protons is not None and self.charge is not None and self.number_of_magnesium_ions is not None:
            return '%s [H%d Z%d Mg%d]' % (self.name, self.number_of_protons or 0, self.charge or 0, self.number_of_magnesium_ions or 0)
        else:
            return '%s' % self.name
    
    def __eq__(self, other):
        """Enable == checking.
        
        Only checks name, protons, charge, and nMg.
        """
        return (str(self.name) == str(other.name) and
                self.number_of_protons == other.number_of_protons and
                self.charge == other.charge and
                self.number_of_magnesium_ions == other.number_of_magnesium_ions)
    
    def __hash__(self):
        """We are HASHABLE!
        
        Note that the hash depends on the same attributes that are checked for equality.
        """
        return hash((self.name, self.number_of_protons, self.charge, self.number_of_magnesium_ions))


class GroupsData(object):
    """Contains data about all groups."""
    
    ORIGIN = Group('Origin', 'Origin', number_of_hydrogens=0, charge=0, number_of_magnesium_ions=0)
    
    # Phosphate groups need special treatment, so they are defined in code...
    # TODO(flamholz): Define them in the groups file.
    
    # each tuple contains: (name, description, nH, charge, nMg, is_default)
    
    phosphate_groups = [('initial H0', '-OPO3-', 0, -1, 0, True),
                        ('initial H1', '-OPO3-', 1, 0, 0, False),
                        ('middle H0', '-OPO2-', 0, -1, 0, True),
                        ('middle H1', '-OPO2-', 1, 0, 0, False),
                        ('final H0', '-OPO3', 0, -2, 0, True),
                        ('final H1', '-OPO3', 1, -1, 0, False),
                        ('final H2', '-OPO3', 2,  0, 0, False),
                        ('initial chain H0', '-OPO3-OPO2-', 0, -2, 0, True),
                        ('initial chain H1', '-OPO3-OPO2-', 1, -1, 0, False),
                        ('initial chain H2', '-OPO3-OPO2-', 2, 0, 0, False),
                        ('initial chain Mg1', '-OPO3-OPO2-', 0, 0, 1, False),
                        ('middle chain H0', '-OPO2-OPO2-', 0, -2, 0, True),
                        ('middle chain H1', '-OPO2-OPO2-', 1, -1, 0, False),
                        ('middle chain H2', '-OPO2-OPO2-', 2, 0, 0, False),
                        ('middle chain Mg1', '-OPO2-OPO2-', 0, 0, 1, False),
                        ('ring initial H0', 'ring -OPO3-', 0, -1, 0, True),
                        ('ring initial H1', 'ring -OPO3-', 1, 0, 0, False),
                        ('ring initial chain H0', 'ring -OPO3-OPO2-', 0, -2, 0, True),
                        ('ring initial chain H1', 'ring -OPO3-OPO2-', 1, -1, 0, False),
                        ('ring initial chain H2', 'ring -OPO3-OPO2-', 2, 0, 0, False),
                        ('ring middle chain H0', 'ring -OPO2-OPO2-', 0, -2, 0, True),
                        ('ring middle chain H1', 'ring -OPO2-OPO2-', 1, -1, 0, False),
                        ('ring middle chain H2', 'ring -OPO2-OPO2-', 2, 0, 0, False),
                        ('ring initial chain Mg1', 'ring -OPO2-OPO2-', 0, 0, 1, False)]
    
    PHOSPHATE_GROUPS = []
    PHOSPHATE_DICT = {}
    DEFAULTS = {}
    for name, desc, number_of_protons, charge, number_of_magnesium_ions, is_default in phosphate_groups:
        group = Group(name, desc, number_of_protons, charge, number_of_magnesium_ions)
        PHOSPHATE_GROUPS.append(group)
        PHOSPHATE_DICT[name] = group
        if is_default:
            DEFAULTS[desc] = group

    RING_PHOSPHATES_TO_MGS = ((PHOSPHATE_DICT['ring initial chain H0'], PHOSPHATE_DICT['ring initial chain Mg1']),)    
    MIDDLE_PHOSPHATES_TO_MGS = ((PHOSPHATE_DICT['initial chain H0'], PHOSPHATE_DICT['initial chain Mg1']),)    
    FINAL_PHOSPHATES_TO_MGS = ((PHOSPHATE_DICT['middle chain H0'], PHOSPHATE_DICT['middle chain Mg1']),)
    
    def __init__(self, groups, transformed=False):
        """Construct GroupsData.
        
        Args:
            groups: a list of Group objects.
        """
        self.transformed = transformed
        self.groups = groups
        self.all_groups = self._filter_groups(self.groups)
        self.all_group_names = [str(g) for g in self.all_groups]
        self.all_group_hydrogens = np.array([g.number_of_protons or 0 for g in self.all_groups])
        self.all_group_charges = np.array([g.charge or 0 for g in self.all_groups])
        self.all_group_mgs = np.array([g.number_of_magnesium_ions or 0 for g in self.all_groups])

        if self.transformed:
            # find the unique group names (ignoring nH, z, nMg)
            # note that Group.name is does not contain these values,
            # unlike Group.__str__() which concatenates the name and the nH, z, nMg
            self.biochemical_group_names = []
            for group in self.all_groups:
                if group.name not in self.biochemical_group_names:
                    self.biochemical_group_names.append(group.name)

    @property
    def count(self):
        return len(self.all_groups)

    @classmethod
    def _filter_groups(cls, groups):
        all_groups = []
        
        for group in groups:
            # Expand phosphate groups.
            if group.is_phosphate():
                all_groups.extend(cls.PHOSPHATE_GROUPS)
            else:
                all_groups.append(group)
        
        # Add the origin.
        all_groups.append(GroupsData.ORIGIN)
        return all_groups
    
    @staticmethod
    def _convert_focal_atoms(focal_atoms_str):
        if not focal_atoms_str:
            return _AllAtomsSet()
        if focal_atoms_str.lower().strip() == 'none':
            return set()
        
        return set([int(c) for c in focal_atoms_str.split('|')])
    
    @staticmethod
    def from_file(filename, transformed=False):
        """Factory that initializes a GroupData from a CSV file."""
        assert filename
        if isinstance(filename, str):
            logging.debug('Reading the list of groups from %s ... ' % filename)
            fp = open(filename, 'r')
        else:
            fp = filename
        list_of_groups = []
        
        gid = 0
        for line_num, row in enumerate(csv.DictReader(fp)):
            logging.debug('Reading group definition for ' + row['NAME'])
            if row.get('SKIP', False):
                logging.debug('Skipping group %s', row.get('NAME'))
                continue
            
            try:
                group_name = row['NAME']
                protons = int(row['PROTONS'])
                charge = int(row['CHARGE'])
                mgs = int(row['MAGNESIUMS'])
                smarts = row['SMARTS']
                focal_atoms = FocalSet(row['FOCAL_ATOMS'])
                #remark = row['REMARK']
                
                # Check that the smarts are good.
                if not Molecule.verify_smarts(smarts):
                    raise GroupsDataError('Cannot parse SMARTS from line %d: %s' % (line_num, smarts))
                
                group = Group(gid, group_name, protons, charge, mgs, str(smarts), focal_atoms)
                list_of_groups.append(group)
            except KeyError as e:
                logging.error(e)
                raise GroupsDataError('Failed to parse row.')
            except ValueError as e:
                logging.error(e)
                raise GroupsDataError('Wrong number of columns (%d) in one of the rows in %s: %s' %
                                      (len(row), filename, str(row)))
            
            gid += 1
        logging.debug('Done reading groups data.')
        
        return GroupsData(list_of_groups, transformed)    

    @staticmethod
    def from_database(db, filename=None, transformed=False):
        """Factory that initializes a GroupData from a DB connection.
        
        Args:
            db: a Database object.
            filename: an optional filename to load data from when
                it's not in the DB. Will write to DB if reading from file.
        
        Returns:
            An initialized GroupsData object.
        """
        logging.debug('Reading the list of groups from the database.')
        
        if not db.DoesTableExist('groups'):
            if filename:
                groups_data = GroupsData.from_file(filename)
                groups_data.to_database(db)
                return groups_data
            else:
                raise Exception('Cannot initialize GroupsData, no file was '
                                'provided and the database does not contain '
                                'the information either')
        
        # Table should exist.
        list_of_groups = []
        for row in db.Execute('SELECT * FROM groups'):
            (gid, group_name, protons, charge, nMg, smarts, focal_atom_set, unused_remark) = row
            try:
                focal_atoms = FocalSet(focal_atom_set)
            except ValueError as e:
                raise ValueError('Group #%d (%s): %s' % (gid, group_name, str(e)))
            list_of_groups.append(Group(gid, group_name, protons, charge, nMg, str(smarts), focal_atoms))
        logging.debug('Done reading groups data.')
        
        return GroupsData(list_of_groups, transformed)
    
    def to_database(self, db):
        """Write the GroupsData to the database."""
        logging.debug('Writing GroupsData to the database.')
        
        db.CreateTable('groups', 'gid INT, name TEXT, protons INT, charge INT, nMg INT, smarts TEXT, focal_atoms TEXT, remark TEXT')
        for group in self.groups:
            focal_atom_str = str(group.focal_atoms)
            db.Insert('groups', [group.id, group.name, group.number_of_protons, group.charge,
                                 group.nMg, group.smarts, focal_atom_str, ''])

        logging.debug('Done writing groups data into database.')

    def index(self, gr):
        try:
            return self.all_groups.index(gr)
        except ValueError:
            raise ValueError('group %s is not defined' % str(gr))

    @property
    def group_names(self):
        if self.transformed:
            return self.biochemical_group_names
        else:
            return self.all_group_names


class GroupDecomposition(object):
    """Class representing the group decomposition of a molecule."""
    
    def __init__(self, groups_data, mol, groups, unassigned_nodes):
        self.groups_data = groups_data
        self.mol = mol
        self.groups = groups
        self.unassigned_nodes = unassigned_nodes
    
    def ToTableString(self):
        """Returns the decomposition as a tabular string."""
        spacer = '-' * 50 + '\n'
        l = ['%30s | %2s | %2s | %3s | %s\n' % ("group name", "nH", "z", "nMg", "nodes"),
             spacer]
                
        for group, node_sets in self.groups:
            if group.number_of_protons is None and group.charge is None and group.nMg is None:
                for n_set in node_sets:
                    s = '%30s |    |    |     | %s\n' % \
                        (group.name, ','.join([str(i) for i in n_set]))
                    l.append(s)
            else:
                for n_set in node_sets:
                    s = '%30s | %2d | %2d | %2d | %s\n' % \
                        (group.name, group.number_of_protons or 0, group.charge or 0, group.nMg or 0,
                         ','.join([str(i) for i in n_set]))
                    l.append(s)

        if self.unassigned_nodes:
            l.append('\nUnassigned nodes: \n')
            l.append('%10s | %3s | %2s | %10s | %10s\n' %
                     ('index', 'an', 'el', 'valence', 'charge'))
            l.append(spacer)
            
            all_atoms = self.mol.GetAtoms()
            for i in self.unassigned_nodes:
                a = all_atoms[i]
                l.append('%10d | %3d | %2s | %10d | %10d\n' %
                         (i, a.GetAtomicNum(), Molecule.GetSymbol(a.GetAtomicNum()),
                          a.GetHvyValence(), a.GetFormalCharge()))
        return ''.join(l)

    def __str__(self):
        """Convert the groups to a string."""        
        group_strs = []
        for group, node_sets in self.non_empty_groups():
            if group.number_of_protons is None and group.charge is None and group.nMg is None:
                group_strs.append('%s x %d' % (group.name, len(node_sets)))
            else:
                group_strs.append('%s [H%d %d %d] x %d' % 
                    (group.name, group.number_of_protons, group.charge, group.nMg,
                     len(node_sets)))
        return " | ".join(group_strs)
    
    def __len__(self):
        counter = 0
        for _group, node_sets in self.non_empty_groups():
            counter += len(node_sets)
        return counter
    
    def as_vector(self):
        """Return the group in vector format.
        
        Note: self.groups contains an entry for *all possible* groups, which is
        why this function returns consistent values for all compounds.
        """
        group_vec = GroupVector(self.groups_data)
        for i, (unused_group, node_sets) in enumerate(self.groups):
            group_vec[i] = len(node_sets)
        group_vec[-1] = 1 # The origin
        return group_vec

    @property
    def non_empty_groups(self):
        """Generator for non-empty groups."""
        for group, node_sets in self.groups:
            if node_sets:
                yield group, node_sets

    @property
    def unassigned_atoms(self):
        """Generator for unassigned atoms."""
        for i in self.unassigned_nodes:
            yield self.mol.GetAtoms()[i], i

    @property
    def sparse_representation(self):
        """Returns a dictionary representation of the group.
        
        TODO(flamholz): make this return some custom object.
        """
        return dict((group, node_sets) for group, node_sets in self.non_empty_groups)

    @property
    def net_charge(self):
        """Returns the net charge."""
        return self.as_vector().net_charge

    @property
    def number_of_protons(self):
        """Returns the number of hydrogens."""
        return self.as_vector().number_of_protons

    @property
    def number_of_magnesium_ions(self):
        """Returns the number of Mg2+ ions."""
        return self.as_vector().number_of_magnesium_ions

    @property
    def count_groups(self):
        """Returns the total number of groups in the decomposition."""
        return sum([len(gdata[-1]) for gdata in self.groups])

    @property
    def pseudoisomer_vectors(self):
        
        def distribute(total, num_slots):
            """
                Returns:
                    a list with all the distinct options of distributing 'total' balls
                    in 'num_slots' slots.
                
                Example:
                    distribute(3, 2) = [[0, 3], [1, 2], [2, 1], [3, 0]]
            """
            if num_slots == 1:
                return [[total]]
            
            if total == 0:
                return [[0] * num_slots]
            
            all_options = []
            for i in range(total+1):
                for opt in distribute(total-i, num_slots-1):
                    all_options.append([i] + opt)
                    
            return all_options
        
        def multi_distribute(total_slots_pairs):
            """
                Returns:
                    similar to distribute, but with more constraints on the sub-totals
                    in each group of slots. Every pair in the input list represents
                    the subtotal of the number of balls and the number of available balls for them.
                    The total of the numbers in these slots will be equal to the subtotal.
                
                Example:
                    multi_distribute([(1, 2), (2, 2)]) =
                    [[0, 1, 0, 2], [0, 1, 1, 1], [0, 1, 2, 0], [1, 0, 0, 2], [1, 0, 1, 1], [1, 0, 2, 0]]
                    
                    in words, the subtotal of the two first slots must be 1, and the subtotal
                    of the two last slots must be 2.
            """
            multilist_of_options = []
            for (total, num_slots) in total_slots_pairs:
                multilist_of_options.append(distribute(total, num_slots))
        
            return [sum(x) for x in itertools.product(*multilist_of_options)]
        
        """Returns a list of group vectors, one per pseudo-isomer."""    
        if not self.count_groups:
            logging.debug('No groups in this decomposition, not calculating pseudoisomers.')
            return []
        
        # A map from each group name to its indices in the group vector.
        # Note that some groups appear more than once (since they can have
        # multiple protonation levels).
        group_name_to_index = {}

        # 'group_name_to_count' is a map from each group name to its number of appearances in 'mol'
        group_name_to_count = {}
        for i, gdata in enumerate(self.groups):
            group, node_sets = gdata
            group_name_to_index.setdefault(group.name, []).append(i)
            group_name_to_count[group.name] = group_name_to_count.get(group.name, 0) + len(node_sets)
        
        index_vector = [] # maps the new indices to the original ones that are used in groupvec

        # A list of per-group pairs (count, # possible protonation levels).
        total_slots_pairs = [] 

        for group_name, groupvec_indices in group_name_to_index.items():
            index_vector += groupvec_indices
            total_slots_pairs.append((group_name_to_count[group_name],
                                      len(groupvec_indices)))

        # generate all possible assignments of protonations. Each group can appear several times, and we
        # can assign a different protonation level to each of the instances.
        groupvec_list = []
        for assignment in multi_distribute(total_slots_pairs):
            v = [0] * len(index_vector)
            for i in range(len(v)):
                v[index_vector[i]] = assignment[i]
            v += [1]  # add 1 for the 'origin' group
            groupvec_list.append(GroupVector(self.groups_data, v))
        return groupvec_list


class GroupDecomposer(object):
    """Decomposes compounds into their constituent groups."""
    
    def __init__(self, groups_data):
        """Construct a GroupDecomposer.
        
        Args:
            groups_data: a GroupsData object.
        """
        self.groups_data = groups_data

    @staticmethod
    def from_file(filename):
        """Factory that initializes a GroupDecomposer from a CSV file."""
        assert filename
        gd = GroupsData.from_file(filename)
        return GroupDecomposer(gd)
    
    @staticmethod
    def from_database(db, filename=None):
        """Factory that initializes a GroupDecomposer from the database.
        
        Args:
            db: a Database object.
            filename: an optional filename to load data from when
                it's not in the DB. Will write to DB if reading from file.
        
        Returns:
            An initialized GroupsData object.
        """
        assert db
        gd = GroupsData.from_database(db, filename)
        return GroupDecomposer(gd)

    @staticmethod
    def _RingedPChainSmarts(length):
        return ''.join(['[C,S][O;R1]', '[P;R1](=O)([OH,O-])[O;R1]' * length, '[C,S]'])

    @staticmethod
    def _InternalPChainSmarts(length):
        return ''.join(['[C,S][O;R0]', '[P;R0](=O)([OH,O-])[O;R0]' * length, '[C,S]'])
    
    @staticmethod
    def _TerminalPChainSmarts(length):
        return ''.join(['[OH,O-]', 'P(=O)([OH,O-])O' * length, '[C,S]'])

    @staticmethod
    def attach_mg_to_phosphate_chain(mol, chain_map, assigned_mgs):
        """Attaches Mg2+ ions the appropriate groups in the chain.
        
        Args:
            mol: the molecule.
            chain_map: the groups in the chain.
            assigned_mgs: the set of Mg2+ ions that are already assigned.
        
        Returns:
            The updated list of assigned Mg2+ ions. 
        """
        # For each Mg2+ we see, we attribute it to a phosphate group if
        # possible. We prefer to assign it to a terminal phosphate, but otherwise 
        # we assign it to a 'middle' group when there are 2 of them.
        def add_mg(p_group, pmg_group, mg):
            node_set = chain_map[p_group].pop(0)
            mg_index = mg[0]
            node_set.add(mg_index)
            assigned_mgs.add(mg_index)
            chain_map[pmg_group].append(node_set)
        
        all_pmg_groups = (GroupsData.FINAL_PHOSPHATES_TO_MGS +
                          GroupsData.MIDDLE_PHOSPHATES_TO_MGS +
                          GroupsData.RING_PHOSPHATES_TO_MGS)
        for mg in mol.find_smarts('[Mg+2]'):
            if mg[0] in assigned_mgs:
                continue
            
            for p_group, pmg_group in all_pmg_groups:
                if chain_map[p_group]:
                    add_mg(p_group, pmg_group, mg)
                    break

        return assigned_mgs

    @staticmethod
    def update_group_map_from_chain(group_map, chain_map):
        """Updates the group_map by adding the chain."""
        for group, node_sets in chain_map.items():
            group_map.get(group, []).extend(node_sets)
        return group_map

    @staticmethod
    def find_phosphate_chains(mol, max_length=4, ignore_protonations=False):
        """
        Chain end should be 'OC' for chains that do not really end, but link to carbons.
        Chain end should be '[O-1,OH]' for chains that end in an hydroxyl.
    
        Args:
            mol: the molecule to decompose.
            max_length: the maximum length of a phosphate chain to consider.
            ignore_protonations: whether or not to ignore protonation values.
        
        Returns:
            A list of 2-tuples (phosphate group, # occurrences).
        """
        group_map = dict((pg, []) for pg in GroupsData.PHOSPHATE_GROUPS)
        v_charge = [a.GetFormalCharge() for a in mol.atoms]
        assigned_mgs = set()
        
        def pop_phosphate(pchain, p_size):
            if len(pchain) < p_size:
                raise Exception('trying to pop more atoms than are left in the pchain')
            phosphate = pchain[0:p_size]
            charge = sum(v_charge[i] for i in phosphate)
            del pchain[0:p_size]
            return set(phosphate), charge
            
        def add_group(chain_map, group_name, charge, atoms):
            default = GroupsData.DEFAULTS[group_name]
            
            if ignore_protonations:
                chain_map[default].append(atoms)
            else:
                # NOTE(flamholz): We rely on the default number of magnesiums being 0 (which it is).
                hydrogens = default.number_of_protons + charge - default.charge
                group = Group(default.id, group_name, hydrogens, charge, default.number_of_magnesium_ions)
                if group not in chain_map:
                    logger.warning('This protonation (%d) level is not allowed for terminal phosphate groups.' %
                                   hydrogens)
                    logger.warning('Using the default protonation level (%d) for this name ("%s").' %
                                   (default.hydrogens, default.name))
                    raise GroupDecompositionError('The group %s cannot have nH = %d' % (group_name, hydrogens))
                else:
                    chain_map[group].append(atoms)
        
        # For each allowed length
        for length in range(1, max_length + 1):
            # Find internal phosphate chains (ones in the middle of the molecule).
            smarts_str = GroupDecomposer._RingedPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.items())
            for pchain in mol.find_smarts(smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop()  # Lose the last carbon
                working_pchain.pop(0)  # Lose the first carbon
                
                if length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 5)
                    add_group(chain_map, 'ring -OPO3-', charge, atoms)                    
                else:
                    atoms, charge = pop_phosphate(working_pchain, 9)
                    add_group(chain_map, 'ring -OPO3-OPO2-', charge, atoms)
                
                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, 'ring -OPO2-OPO2-', charge, atoms)
            
            assigned_mgs = GroupDecomposer.attach_mg_to_phosphate_chain(mol, chain_map, assigned_mgs)
            GroupDecomposer.update_group_map_from_chain(group_map, chain_map)

            # Find internal phosphate chains (ones in the middle of the molecule).
            smarts_str = GroupDecomposer._InternalPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.items())
            for pchain in mol.find_smarts(smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop() # Lose the last carbon
                working_pchain.pop(0) # Lose the first carbon
                
                if length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 5)
                    add_group(chain_map, '-OPO3-', charge, atoms)                    
                else:
                    atoms, charge = pop_phosphate(working_pchain, 9)
                    add_group(chain_map, '-OPO3-OPO2-', charge, atoms)
                
                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, '-OPO2-OPO2-', charge, atoms)
            
            assigned_mgs = GroupDecomposer.attach_mg_to_phosphate_chain(mol, chain_map, assigned_mgs)
            GroupDecomposer.update_group_map_from_chain(group_map, chain_map)
            
            # Find terminal phosphate chains.
            smarts_str = GroupDecomposer._TerminalPChainSmarts(length)
            chain_map = dict((k, []) for (k, _) in group_map.items())
            for pchain in mol.find_smarts(smarts_str):
                working_pchain = list(pchain)
                working_pchain.pop()  # Lose the carbon
                
                atoms, charge = pop_phosphate(working_pchain, 5)
                add_group(chain_map, '-OPO3', charge, atoms)
                
                if not length % 2:
                    atoms, charge = pop_phosphate(working_pchain, 4)
                    add_group(chain_map, '-OPO2-', charge, atoms)
                
                while working_pchain:
                    atoms, charge = pop_phosphate(working_pchain, 8)
                    add_group(chain_map, '-OPO2-OPO2-', charge, atoms)
                
            assigned_mgs = GroupDecomposer.attach_mg_to_phosphate_chain(mol, chain_map, assigned_mgs)
            GroupDecomposer.update_group_map_from_chain(group_map, chain_map)

        return [(pg, group_map[pg]) for pg in GroupsData.PHOSPHATE_GROUPS]

    def create_empty_group_decomposition(self):
        emptymol = Molecule.from_smiles("")
        decomposition = self.decompose(emptymol, ignore_protonations=True, strict=False)
        for i, (group, _node_sets) in enumerate(decomposition.groups):
            decomposition.groups[i] = (group, [])
        return decomposition

    def decompose(self, mol, ignore_protonations=False, strict=False):
        """
        Decompose a molecule into groups.
        
        The flag 'ignore_protonations' should be used when decomposing a compound with lacing protonation
        representation (for example, the KEGG database doesn't posses this information). If this flag is
        set to True, it overrides the '(C)harge sensitive' flag in the groups file (i.e. - *PC)
        
        Args:
            mol: the molecule to decompose.
            ignore_protonations: whether to ignore protonation levels.
            strict: whether to assert that there are no unassigned atoms.
        
        Returns:
            A GroupDecomposition object containing the decomposition.
        """
        unassigned_nodes = set(range(len(mol)))
        groups = []
        
        def _add_correction(group, count):
            l = [set() for _ in range(count)]
            groups.append((group, l))
        
        for group in self.groups_data.groups:
            # Phosphate chains require a special treatment
            if group.is_phosphate():
                pchain_groups = None
                if group.ignore_charges() or ignore_protonations:
                    pchain_groups = self.find_phosphate_chains(mol, ignore_protonations=True)
                elif group.charge_sensitive():
                    pchain_groups = self.find_phosphate_chains(mol, ignore_protonations=False)
                else:
                    raise MalformedGroupDefinitionError(
                        'Unrecognized phosphate wildcard: %s' % group.name)
                
                for phosphate_group, group_nodesets in pchain_groups:
                    current_groups = []
                    
                    for focal_set in group_nodesets:
                        if focal_set.issubset(unassigned_nodes):
                            # Check that the focal-set doesn't override an assigned node
                            current_groups.append(focal_set)
                            unassigned_nodes = unassigned_nodes - focal_set
                    groups.append((phosphate_group, current_groups))
            elif group.is_coded_correction():
                _add_correction(group, group.correction(mol))
            # Not a phosphate group or expanded correction.
            else:
                # TODO: if the 'ignore_protonation' flag is True, this should always
                # use the pseudogroup with the lowest nH in each category regardless
                # of the hydrogens in the given Mol.
                current_groups = []
                for nodes in mol.find_smarts(group.smarts):
                    try:
                        focal_set = group.focal_set(nodes)
                    except IndexError:
                        logging.error('Focal set for group %s is out of range: %s'
                                      % (str(group), str(group.focal_atoms)))
                        sys.exit(-1)

                    # check that the focal-set doesn't override an assigned node
                    if focal_set.issubset(unassigned_nodes): 
                        current_groups.append(focal_set)
                        unassigned_nodes = unassigned_nodes - focal_set
                groups.append((group, current_groups))
        
        # Ignore the hydrogen atoms when checking which atom is unassigned
        for nodes in mol.find_smarts('[H]'):
            unassigned_nodes = unassigned_nodes - set(nodes)
        
        decomposition = GroupDecomposition(self.groups_data, mol, groups, unassigned_nodes)
        
        if strict and decomposition.unassigned_nodes:
            raise GroupDecompositionError('Unable to decompose %s into groups.' % mol.title, decomposition)
        
        return decomposition


def inchi_to_group_vector(group_decomposer, inchi):
    try:
        mol = Molecule.from_inchi(str(inchi))
    except OpenBabelError:
        raise GroupDecompositionError('cannot convert InChI to Molecule')

    # mol.RemoveHydrogens()
    decomposition = group_decomposer.decompose(mol, ignore_protonations=False, strict=True)

    # nH = decomposition.Hydrogens()
    # charge = decomposition.NetCharge()
    # nMg = decomposition.Magnesiums()
    return decomposition.as_vector()


def smiles_to_group_vector(group_decomposer, smiles):
    try:
        mol = Molecule.from_smiles(str(smiles))
    except OpenBabelError:
        raise GroupDecompositionError('cannot convert InChI to Molecule')

    # mol.RemoveHydrogens()
    decomposition = group_decomposer.decompose(mol, ignore_protonations=False, strict=True)

    # nH = decomposition.Hydrogens()
    # charge = decomposition.NetCharge()
    # nMg = decomposition.Magnesiums()
    return decomposition.as_vector()


def init_groups_data():
    with open(GROUP_CSV, 'r') as groups_file:
        return GroupsData.from_file(groups_file, transformed=False)


