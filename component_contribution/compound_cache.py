import os
import json
import gzip
import numpy
import logging
from pandas import DataFrame, read_csv

from component_contribution.core.compound import Compound

from component_contribution.core.singleton import Singleton


logger = logging.getLogger(__name__)

base_path = os.path.split(os.path.realpath(__file__))[0]

# Input Files:
# original version of the KEGG compound file 
OLD_COMPOUND_JSON_FNAME = os.path.join(base_path, '../data/equilibrator_compounds.json.gz')

# a CSV file with additional names and InChIs (mostly compounds missing from KEGG and added manually)
KEGG_ADDITIONS_TSV_FNAME = os.path.join(base_path, '../data/kegg_additions.tsv')

# Files created by this module:
# names and InChIs only
KEGG_COMPOUND_JSON_FNAME = os.path.join(base_path, '../data/kegg_compounds.json.gz') 

# names, InChIs and pKa data
DEFAULT_CACHE_FNAME = os.path.join(base_path, '../cache/compounds.json.gz') 


class CompoundCache(Singleton):
    """
    CompoundCache is a singleton that handles caching of Compound objects for the component-contribution package.

    The Compounds are retrieved by their ID (which is the KEGG ID in most cases). The first time a Compound is
    requested, it is obtained from the relevant database and a Compound object is created (this takes a while because
    it usually involves internet communication and then invoking the ChemAxon plugin for calculating the pKa values
    for that structure).

    Any further request for the same Compound ID will draw the object from the cache. When the method dump() is called,
    all cached data is written to a file that will be loaded in future python sessions.

    Attributes
    ----------
    cache_file_name : str
        The name of the file where the cache is dumped (default '.compound.cache').
    _data : DataFrame
        The internal cache storage.
    need_to_update_cache_file : bool
        Whether the cache needs update.
    inchi_key2compound_id : dict
        Mapping from InChi Keys to compound ids

    """

    COLUMNS = ["inchi", "atom_bag", "compound_object"]

    def __init__(self, cache_file_name=None):
        if cache_file_name is None:
            cache_file_name = ".compound.cache"
        self.cache_file_name = cache_file_name
        self._data = DataFrame(columns=self.COLUMNS)
        self.need_to_update_cache_file = False
        self.inchi_key2compound_id = {}
        self.load()

    def __contains__(self, item):
        if isinstance(item, str):
            return item in self._data.index or item in self.inchi_key2compound_id
        else:
            return False

    @property
    def compound_ids(self):
        return self._data.index.tolist()
    
    def load(self):
        if os.path.exists(self.cache_file_name):
            data = read_csv(self.cache_file_name)
            for compound_id, row in data.iterrows():
                compound = Compound.from_json_dict(json.loads(row.compound_object))
                self._data.loc[compound_id] = [row.inchi, compound.atom_bag, compound]
                if compound.inchi:
                    self.inchi_key2compound_id[compound.inchi_key] = compound_id

    def dump(self):
        if self.need_to_update_cache_file:
            to_dump = DataFrame(columns=["inchi", "compound_object"])
            for compound_id, row in self._data.iterrows():
                compound_json = json.dumps(row.compound_object.to_json_dict())
                to_dump.loc[compound_id] = [row.inchi, compound_json]

            to_dump.to_csv(self.cache_file_name)
            self.need_to_update_cache_file = False
        
    def get_compound(self, compound_id):
        if compound_id in self.inchi_key2compound_id:
            compound_id = self.inchi_key2compound_id[compound_id]

        if compound_id not in self._data.index:
            logger.debug('Cache missing: %s' % str(compound_id))
            comp = Compound.from_database(compound_id)
            self.add(comp)

        logger.debug('Cache hit: %s' % str(compound_id))
        return self._data.loc[compound_id, 'compound_object']

    def remove(self, compound_id):
        if compound_id in self._data.index:
            compound = self._data.loc[compound_id, 'compound_object']
            if compound.inchi_key in self.inchi_key2compound_id:
                del self.inchi_key2compound_id[compound.inchi_key]
            self._data.drop(compound_id)
        else:
            raise KeyError(compound_id)
    
    def add(self, compound: Compound):
        if compound.compound_id in self._data.index:
            raise KeyError(compound.compound_id)
        if compound.inchi_key in self.inchi_key2compound_id:
            raise KeyError("%s is the same as %s" %
                           (compound.compound_id, self.inchi_key2compound_id[compound.inchi_key.compound.inchi_key]))

        self._data.loc[compound.compound_id] = [compound.inchi, compound.atom_bag, compound]
        self.inchi_key2compound_id[compound.inchi_key] = compound.compound_id
        self.need_to_update_cache_file = True
            
    def get_element_matrix(self, compound_ids):
        if type(compound_ids) == str:
            compound_ids = [compound_ids]
        # gather the "atom bags" of all compounds in a list 'atom_bag_list'
        elements = set()
        atom_bag_list = []
        for compound_id in compound_ids:
            comp = self.get_compound(compound_id)
            atom_bag = comp.atom_bag
            if atom_bag is not None:
                elements = elements.union(atom_bag.keys())
            atom_bag_list.append(atom_bag)
        elements.discard('H')  # don't balance H (it's enough to balance e-)
        elements = sorted(elements)
        
        # create the elemental matrix, where each row is a compound and each column is an element (or e-)
        element_matrix = numpy.matrix(numpy.zeros((len(atom_bag_list), len(elements))))
        for i, atom_bag in enumerate(atom_bag_list):
            if atom_bag is None:
                element_matrix[i, :] = numpy.nan
            else:
                for j, elem in enumerate(elements):
                    element_matrix[i, j] = atom_bag.get(elem, 0)
        return elements, element_matrix


def load_kegg_compounds_from_json(file_name, cache: CompoundCache):
    compounds = json.load(gzip.open(file_name, 'r'))
    for compound in compounds:
        cache.add(Compound.from_json_dict(compound))
