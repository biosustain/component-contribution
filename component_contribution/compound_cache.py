# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import atexit
import json
import logging
import sqlite3
from collections import defaultdict
from contextlib import ExitStack

import pandas as pd
import pandas.io.sql as pd_sql
from importlib_resources import path

import component_contribution.cache
from .compound import Compound, MicroSpecie


logger = logging.getLogger(__name__)
file_manager = ExitStack()
atexit.register(file_manager.close)

DEFAULT_CACHE_FNAME = str(file_manager.enter_context(
    path(component_contribution.cache, "compounds.sqlite")))


class CompoundCache(object):
    """
    CompoundCache is a singleton that handles caching of Compound objects for
    the component-contribution package.  The Compounds are retrieved by their
    ID (e.g., KEGG COMPOUND ID, ChEBI Id or HMDB in most cases) or InChI Key.
    The first time a Compound is requested, it is obtained from the relevant
    database and a Compound object is created (this takes a while because it
    usually involves internet communication and then invoking the ChemAxon
    plugin for calculating the pKa values for that structure). Any further
    request for the same Compound ID will draw the object from the cache. When
    the method dump() is called, all cached data is written to a file that will
    be loaded in future python sessions.
    """

    COLUMNS = [
        "inchi", "name", "atom_bag", "p_kas", "smiles", "major_ms",
        "number_of_protons", "charges"]

    def __init__(self, cache_file_name=DEFAULT_CACHE_FNAME):
        # connect to the SQL database
        self.con = sqlite3.connect(cache_file_name)

        # a lookup table for compound IDs
        self._compound_id_to_inchi_key = dict()

        # a lookup table for InChIKeys
        self._inchi_key_to_compound_ids = defaultdict(set)

        # a lookup table to microspecies (charge, nH, energy difference)
        self._inchi_key_to_species = dict()

        # an internal cache for Compound objects that have already been created
        self.compound_dict = dict()

        # a flag telling the Cache that it needs to rewrite itself in the
        # filesystem. updated any time a new compound is cached.
        self.requires_update = False

        if pd_sql.has_table('compounds', self.con):
            compound_df = pd.read_sql('SELECT * FROM compounds', self.con)
            specie_df = pd.read_sql('SELECT * FROM microspecies', self.con)
            self.from_data_frames(compound_df, specie_df)
        elif pd_sql.has_table('compound_cache', self.con):
            df = pd.read_sql('SELECT * FROM compound_cache', self.con)
            self.from_data_frame(df)
            self.requires_update = True
            self.dump()
        else:
            self._data = pd.DataFrame(columns=self.COLUMNS)
            self._data.index.name = 'inchi_key'
            self.requires_update = True
            self.dump()

#    def get_species_frame(self):
#        """
#            Create a pandas.DataFrame that contains all the pseudoisomer info,
#            i.e. the net charge, number of protons, and deltaG0
#        """
#        R_T_log10 = R * default_T * np.log(10)
#        data = []
#        for inchi_key, row in self._data.iterrows():
#            major_ms = row['major_ms']
#            p_kas = row['p_kas']
#            for i, (z, nH) in enumerate(zip(row['charges'],
#                                            row['number_of_protons'])):
#                if i < major_ms:
#                    ddG = sum(p_kas[i:major_ms]) * R_T_log10
#                elif i > major_ms:
#                    ddG = -sum(p_kas[major_ms:i]) * R_T_log10
#                else:
#                    ddG = 0
#
#                data.append((inchi_key, z, nH, 0, ddG))
#
#        species_df = pd.DataFrame(data=data,
#                                  columns=['inchi_key',
#                                           'charge',
#                                           'number_of_protons',
#                                           'number_of_magnesium',
#                                           'relative_dG0_f'])
#        species_df['relative_dG0_f'] = species_df['relative_dG0_f'].round(2)
#        return species_df
#
    def dump(self, cache_file_name=DEFAULT_CACHE_FNAME):
        if self.requires_update:
            compound_df, specie_df = self.to_data_frames()
            pd_sql.to_sql(compound_df, 'compounds', self.con,
                          if_exists='replace')
            pd_sql.to_sql(specie_df, 'microspecies', self.con,
                          if_exists='replace')
            self.requires_update = False
            self.con.commit()

    def _read_cross_refs(self, df):
        for inchi_key, row in df.iterrows():
            for compound_id in row.cross_references.split(";"):
                self._compound_id_to_inchi_key[compound_id] = inchi_key
                self._inchi_key_to_compound_ids[inchi_key].add(compound_id)

    def all_compound_ids(self):
        return sorted(self._compound_id_to_inchi_key.keys())

    def to_data_frames(self):
        """
            Creates DataFrames that can be written to the SQLite database.
            The atom_bag columns needs to be serialized before the export.

            Returns:
                (pandas.DataFrame, pandas.DataFrame)
        """
        compound_df = self._data.copy()

        print(self._inchi_key_to_compound_ids)
        compound_df['cross_references'] = \
            compound_df.index.map(self._inchi_key_to_compound_ids.get)
        compound_df['cross_references'] = \
            compound_df['cross_references'].apply(';'.join)

        compound_df['atom_bag'] = compound_df['atom_bag'].apply(json.dumps)

        specie_df = pd.DataFrame.from_dict(
            [s.to_dict()
             for species in self._inchi_key_to_species.values()
             for s in species])

        return compound_df, specie_df

    def from_data_frame(self, df):
        """
            Legacy code, helping in the transition to several SQL tables
        """
        # read the compound DataFrame and adjust the columns that are not
        # straighforward
        self._data = df.copy()
        self._data.set_index('inchi_key', inplace=True, drop=True)
        self._data['atom_bag'] = self._data['atom_bag'].apply(json.loads)
        self._data['charges'] = self._data['charges'].apply(json.loads)
        self._data['number_of_protons'] = self._data['number_of_protons'].apply(json.loads)
        self._data['p_kas'] = self._data['p_kas'].apply(json.loads)
        self._data['inchi'].fillna('', inplace=True)
        self._data['smiles'].fillna('', inplace=True)

        self._read_cross_refs(self._data)
        self._data.drop('cross_references', axis=1, inplace=True)

        for inchi_key, row in self._data.iterrows():
            self._inchi_key_to_species[inchi_key] = \
                MicroSpecie.from_p_kas(inchi_key, int(row.major_ms),
                                       row.charges,
                                       row.number_of_protons,
                                       row.p_kas)
        self._data.drop(['charges', 'number_of_protons', 'p_kas', 'major_ms'],
                        axis=1, inplace=True)

    def from_data_frames(self, compound_df, specie_df):
        """
            Reads all the cached data from DataFrames

            Arguments:
                compound_df - a compound DataFrame created by to_data_frames()
                specie_df -  a specie DataFrame created by to_data_frames()
        """
        # read the compound DataFrame and adjust the columns that are not
        # straighforward
        self._data = compound_df.copy()
        self._data.set_index('inchi_key', inplace=True, drop=True)
        self._data['atom_bag'] = self._data['atom_bag'].apply(json.loads)
        self._data['inchi'].fillna('', inplace=True)
        self._data['smiles'].fillna('', inplace=True)

        self._read_cross_refs(self._data)
        self._data.drop('cross_references', axis=1, inplace=True)

        # read the microspecie DataFrame into a dictionary
        for inchi_key, group_df in specie_df.groupby('inchi_key'):
            species = []
            for _, row in group_df.sort_values('charge').iterrows():
                ms = MicroSpecie(inchi_key, row.charge, row.number_of_protons,
                                 row.ddG_over_T, bool(row.is_major))
                species.append(ms)
            self._inchi_key_to_species[inchi_key] = species

        for inchi_key in self._data.index:
            if inchi_key not in self._inchi_key_to_species:
                self._inchi_key_to_species[inchi_key] = []

#        for inchi_key, row in self._data.iterrows():
#            self._inchi_key_to_species[inchi_key] = \
#                MicroSpecie.from_p_kas(inchi_key, int(row.major_ms),
#                                       row.charges,
#                                       row.number_of_protons,
#                                       row.p_kas)

    def exists(self, compound_id):
        return ((compound_id in self._compound_id_to_inchi_key) or
                (compound_id in self._data.index))

    def add_cross_link(self, inchi_key, new_compound_id):
        logging.debug('Adding a cross-link from %s to %s' %
                      (new_compound_id, inchi_key))
        cpd = self.get(inchi_key)
        self._compound_id_to_inchi_key[new_compound_id] = cpd.inchi_key
        self.requires_update = True

    def get_compound(self, compound_id, compute_pkas=True):
        # compound exists
        if compound_id in self._compound_id_to_inchi_key:
            logging.debug('Cache hit for %s' % compound_id)
            return self.get(self._compound_id_to_inchi_key[compound_id])
        # compound_id is an InChI Key and exists
        elif compound_id in self._data.index:
            logging.debug('Cache hit for InChiKey %s' % compound_id)
            return self.get(compound_id)
        # compound does not exist.
        else:
            logging.debug('Cache miss, calculating pKas for %s' % compound_id)
            cpd = Compound.get(compound_id, compute_pkas)
            if cpd.inchi_key in self._data.index:
                self.add_cross_link(cpd.inchi_key, compound_id)
            else:
                logging.debug('Adding the new InChiKey to the cache: %s'
                              % cpd.inchi_key)
                self.add(cpd)
            return cpd

    def remove(self, inchi_key):
        if inchi_key not in self._data.index:
            raise KeyError(inchi_key)

        del self.compound_dict[inchi_key]
        del self._inchi_key_to_compound_ids
        self._data.drop(inchi_key)

        for compound_id in self._compound_id_to_inchi_key:
            if self._compound_id_to_inchi_key[compound_id] == inchi_key:
                del self._compound_id_to_inchi_key[compound_id]

    def get(self, inchi_key):
        if inchi_key not in self._data.index:
            raise KeyError(inchi_key)

        if inchi_key in self.compound_dict:
            cpd = self.compound_dict[inchi_key]
        else:
            data = self._data.loc[inchi_key, :]
            species = self._inchi_key_to_species[inchi_key]
            cpd = Compound(
                inchi_key, data.inchi, data.atom_bag, species, data.smiles)
            self.compound_dict[inchi_key] = cpd

        return cpd

    def add(self, cpd):
        if cpd.inchi_key not in self._inchi_key_to_compound_ids:
            self._inchi_key_to_compound_ids[cpd.inchi_key] = set()

        self._inchi_key_to_compound_ids[cpd.inchi_key].add(cpd.compound_id)
        self._compound_id_to_inchi_key[cpd.compound_id] = cpd.inchi_key
        self.compound_dict[cpd.inchi_key] = cpd
        self._inchi_key_to_species[cpd.inchi_key] = cpd.species
        self._data.loc[cpd.inchi_key, :] = [
            cpd.inchi, cpd.name, cpd.atom_bag, cpd.p_kas, cpd.smiles,
            int(cpd.major_microspecies), cpd.number_of_protons, cpd.charges]
        self.requires_update = True

    def get_element_data_frame(self, compound_ids):
        if isinstance(compound_ids, str):
            compound_ids = [compound_ids]

        # gather the "atom bags" of all compounds in a list 'atom_bag_list'
        atom_bags = [self.get_compound(cid).atom_bag or {}
                     for cid in compound_ids]

        elements = sorted(set([e for bag in atom_bags for e in bag.keys()]))

        # create the elemental matrix, where each row is a compound and each
        # column is an element (or e-)
        element_df = pd.DataFrame(index=compound_ids, columns=elements,
                                  dtype=int).fillna(0)

        for cid, atom_bag in zip(compound_ids, atom_bags):
            for elem in elements:
                element_df[elem][cid] = atom_bag.get(elem, 0)

        return element_df


# this is the only place where one should use the constructore.
# we wish to only have one instance of the cache (i.e. use it as a singleton)

# legacy code for transitioning from CSV to SQLite
if False:
    con = sqlite3.connect(DEFAULT_CACHE_FNAME)
    df = pd.read_csv('component_contribution/cache/compounds.csv',
                     index_col=0)
    df.index.name = 'inchi_key'
    df.to_sql('compound_cache', con, if_exists='replace')
    con.commit()
    con.close()

    con = sqlite3.connect(DEFAULT_CACHE_FNAME)
    df = pd.read_sql('SELECT * FROM compound_cache', con)
    df.set_index('inchi_key', inplace=True)
    print(df)

    con.close()

ccache = CompoundCache()


if __name__ == '__main__':
    from .training_data import FullTrainingData
    logger.setLevel(logging.DEBUG)

    # Cache all the compounds that are part of the training data.
    # Calling the constructor here,
    # already invokes queries for all the relevant compounds (since their
    # structure is needed in order to balance the reactions and the pKas
    # are needed in order to perform the reverse Legendre transform).
    td = FullTrainingData()

    ccache.dump()
