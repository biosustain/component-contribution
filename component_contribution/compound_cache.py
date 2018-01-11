import os
import json
import gzip
import numpy
import logging

import time
from pandas import DataFrame, read_csv

from component_contribution.core.compound import Compound

from component_contribution.core.singleton import Singleton


logger = logging.getLogger(__name__)

current_directory = os.path.dirname(__file__)


PERMANENT_CACHE = os.path.join(current_directory, '../data/cached_compounds.json.gz')


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
    compound_id2inchi_key : dict
        Mapping from InChi Keys to compound ids.
    """

    COLUMNS = ["inchi", "atom_bag", "compound_object"]

    def __init__(self, cache_file_name=None, bootstrap=False):
        if cache_file_name is None:
            cache_file_name = ".compound.cache"
        self.cache_file_name = cache_file_name
        self._data = DataFrame(columns=self.COLUMNS)
        self.need_to_update_cache_file = False
        self.compound_id2inchi_key = {}
        self.inchi_key2compound_id = {}
        self.load()
        if bootstrap:
            self.bootstrap()

    def __len__(self):
        return len(self._data)

    def bootstrap(self):
        """
        To be used when the cache does not exist. It loads the existing chemicals from a file.
        """
        with gzip.open(PERMANENT_CACHE, 'rt', encoding='utf-8') as data:
            compounds = json.load(data)
            for compound in compounds:
                compound = Compound.from_json_dict(compound)
                if compound not in self:
                    self.add(compound)

            logger.debug("BOOTSTRAP: added %i compounds to cache" % len(compounds))

        self.dump()

    def __contains__(self, item):
        if isinstance(item, str):
            return item in self._data.index or item in self.compound_id2inchi_key
        if isinstance(item, Compound):
            if item.inchi_key is not None:
                return item.inchi_key in self._data.index
            else:
                return item.compound_id in self._data.index
        else:
            return False

    @property
    def compound_ids(self):
        return self._data.index.tolist()
    
    def load(self):
        def _parse_compound_ids(data):
            if isinstance(data, float):
                return set()
            elif isinstance(data, str):
                return set(data.split(";"))
            else:
                raise TypeError("unexpect type %s" % type(data))

        def _map_compound_index(compound):
            if compound.inchi is not None:
                return compound.inchi_key
            else:
                return compound.compound_id

        if os.path.exists(self.cache_file_name):

            data = read_csv(self.cache_file_name, converters={"compound_ids": _parse_compound_ids})
            new_data = numpy.ndarray((len(data.index), len(self.COLUMNS)), dtype=object)
            existing = dict()

            start = time.time()
            for i, (entry_id, row) in enumerate(data.iterrows()):
                compound = Compound.from_json_dict(json.loads(row.compound_object))

                idx = _map_compound_index(compound)

                if idx in existing:
                    new_data[i] = [numpy.nan, numpy.nan, numpy.nan]
                    continue
                else:
                    existing[idx] = None
                    new_data[i] = [row.inchi, compound.atom_bag, compound]

                if compound.inchi is not None:
                    self.compound_id2inchi_key[compound.compound_id] = compound.inchi_key
                    for cid in row.compound_ids:
                        self.compound_id2inchi_key[cid] = compound.inchi_key
                    self.inchi_key2compound_id[compound.inchi_key] = row.compound_ids

                if i > 0 and i % 1000 == 0:
                    logger.debug("%i entries\t %.2f seconds" % (i, time.time() - start))

            self._data = DataFrame(data=new_data, columns=self.COLUMNS)
            self._data.dropna(inplace=True)
            self._data.index = self._data.compound_object.apply(_map_compound_index)

            logger.debug("Index: %s" % self._data.index.tolist()[:3])

    def dump(self):
        if self.need_to_update_cache_file:
            data = numpy.ndarray((len(self._data.index), 3), dtype=object)

            for i, (compound_id, row) in enumerate(self._data.iterrows()):
                compound_json = json.dumps(row.compound_object.to_json_dict())
                compound_ids = self.inchi_key2compound_id.get(row.compound_object.id, [])

                data[i] = [row.inchi, compound_json, ";".join(compound_ids)]

            to_dump = DataFrame(data, index=self._data.index.tolist(),
                                columns=["inchi", "compound_object", "compound_ids"])

            to_dump.to_csv(self.cache_file_name)
            self.need_to_update_cache_file = False
        
    def get_compound(self, compound_id):
        logger.debug("get_compound(%s)" % compound_id)

        if compound_id is None:
            return Compound.from_database("")

        if compound_id in self.compound_id2inchi_key:
            compound_id = self.compound_id2inchi_key[compound_id]
            compound = self._data.loc[compound_id].compound_object
            logger.debug('found: %s' % compound.id)

        if compound_id not in self._data.index:
            logger.debug('missing: %s' % str(compound_id))
            compound = Compound.from_database(compound_id)
            self.add(compound)
        else:
            compound = self._data.loc[compound_id].compound_object
            logger.debug('found: %s' % compound.id)

        return compound

    def remove(self, compound_id):
        if compound_id in self._data.index:

            compound = self._data.loc[compound_id].compound_object
            if compound.compound_id in self.compound_id2inchi_key:
                del self.compound_id2inchi_key[compound.compound_id]
            if compound.inchi_key in self.inchi_key2compound_id:
                del self.inchi_key2compound_id[compound.inchi_key]

            self._data.drop(compound_id, inplace=True)
        else:
            raise KeyError(compound_id)

    def replace(self, compound):
        if compound.compound_id in self._data.index:
            index = compound.compound_id
        elif compound.inchi_key is not None and compound.inchi_key in self._data.index:
            index = compound.inchi_key
        else:
            raise KeyError("Compound %s is not cached" % compound.id)

        self.remove(index)
        self.add(compound)
    
    def add(self, compound: Compound):
        if compound.inchi_key is None:
            cid = compound.compound_id
        else:
            cid = compound.inchi_key

            if cid in self._data.index:
                self.compound_id2inchi_key[compound.compound_id] = compound.inchi_key
                self.inchi_key2compound_id[cid].add(compound.compound_id)
                return

            if cid in self.inchi_key2compound_id and compound.compound_id in self.inchi_key2compound_id[cid]:
                raise KeyError("%s is the same as %s (%s)" %
                               (compound.compound_id, self.inchi_key2compound_id[cid]), cid in self._data.index)
            else:
                self.inchi_key2compound_id[cid] = {compound.compound_id}
                self.compound_id2inchi_key[compound.compound_id] = cid

        if cid in self._data.index:
            raise KeyError(cid)

        self._data.loc[cid] = [compound.inchi, compound.atom_bag, compound]

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
