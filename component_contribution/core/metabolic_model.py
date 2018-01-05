import copy
import csv
import logging
import re

import numpy as np

from component_contribution.compound_cache import CompoundCache
from component_contribution.core.reaction import Reaction
from component_contribution.exceptions import ParseException

logger = logging.getLogger(__name__)


class MetabolicModel(object):

    def __init__(self, s_matrix, compound_ids, database, rids=None):
        self.s_matrix = s_matrix
        self.s_matrix_mets = copy.copy(s_matrix)
        self.compound_ids = compound_ids
        self.rids = rids
        self.database = database

        assert len(self.compound_ids) == self.s_matrix.shape[0]

        if self.rids is not None:
            assert len(self.rids) == self.s_matrix.shape[1]

        self.ccache = CompoundCache.get_instance()

        # DO THIS FOR OTHER FORMATS AS WELL
        # remove H+ from the stoichiometric matrix if it exists
        while 'C00080' in self.compound_ids:
            i = self.compound_ids.index('C00080')
            self.s_matrix = np.vstack((self.s_matrix[:i, :], self.s_matrix[i + 1:, :]))
            self.compound_ids.pop(i)

        self.dg0 = None
        self.cov_dg0 = None

    def __del__(self):
        if hasattr(self, 'ccache'):
            self.ccache.dump()

    @classmethod
    def from_file(cls, file_name, arrow=None, input_format='kegg', has_reaction_ids=False):
        """
        Reads a file containing reactions in KEGG format
        
        Arguments
        ---------
        file_name : str
            The name of the file to read.
        arrow : str
            The string used as the 'arrow' in each reaction (default: '<=>')
        input_format : str
            The text file format provided ('kegg', 'tsv' or 'csv')
        has_reaction_ids : bool
            A flag indicating if there is a column of reaction IDs (separated from the reaction with whitespaces)
        
        Returns
        -------
        MetabolicModel
        """

        if arrow is None:
            arrow = "<=>"
        with open(file_name, 'r') as input_file:
            if input_format == 'kegg':
                model = MetabolicModel.from_kegg_formulas(input_file.readlines(), arrow, has_reaction_ids)
            elif input_format == "bigg":
                model = MetabolicModel.from_formulas(input_file.readlines(), database='bigg', arrow=arrow,
                                                     has_reaction_ids=has_reaction_ids)
            elif input_format == 'tsv':
                model = MetabolicModel.from_csv(input_file, has_reaction_ids=has_reaction_ids, delimiter='\t')
            elif input_format == 'csv':
                model = MetabolicModel.from_csv(input_file, has_reaction_ids=has_reaction_ids, delimiter=None)
            else:
                raise ValueError("Invalid input format %s (use kegg, tsv or csv)" % input_format)

            return model
    
    @staticmethod
    def from_csv(fd, has_reaction_ids=True, delimiter=None):
        csv_reader = csv.reader(fd, delimiter=delimiter)
        if has_reaction_ids:
            rids = csv_reader.next()
            rids = rids[1:]
        else:
            rids = None
        s_matrix = []
        compound_ids = []
        for i, row in enumerate(csv_reader):
            compound_ids.append(row[0])
            s_matrix.append([float(x) for x in row[1:]])
        s_matrix = np.array(s_matrix)

        return MetabolicModel(s_matrix, compound_ids, rids)

    @staticmethod
    def from_reactions(reactions, format, has_reaction_ids=False):
        if has_reaction_ids:
            reaction_ids = [r.rid for r in reactions]
        else:
            reaction_ids = None

        compound_ids = set()
        for reaction in reactions:
            compound_ids = compound_ids.union(reaction.keys())

        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'compound_ids'.
        compound_ids = sorted(compound_ids)
        s_matrix = np.matrix(np.zeros((len(compound_ids), len(reactions))))
        for i, reaction in enumerate(reactions):
            s_matrix[:, i] = np.matrix(reaction.dense(compound_ids))

        logger.debug('Successfully loaded %d reactions (involving %d unique compounds)' %
                     (s_matrix.shape[1], s_matrix.shape[0]))
        return MetabolicModel(s_matrix, compound_ids, format, reaction_ids)

    @staticmethod
    def from_kegg_reactions(kegg_reactions, has_reaction_ids=False):
        if has_reaction_ids:
            reaction_ids = [r.rid for r in kegg_reactions]
        else:
            reaction_ids = None

        compound_ids = set()
        for reaction in kegg_reactions:
            compound_ids = compound_ids.union(reaction.stoichiometry.keys())
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'compound_ids'.
        compound_ids = sorted(compound_ids)
        s_matrix = np.matrix(np.zeros((len(compound_ids), len(kegg_reactions))))
        for i, reaction in enumerate(kegg_reactions):
            s_matrix[:, i] = np.matrix(reaction.dense(compound_ids))
        
        logger.debug('Successfully loaded %d reactions (involving %d unique compounds)' %
                     (s_matrix.shape[1], s_matrix.shape[0]))
        return MetabolicModel(s_matrix, compound_ids, reaction_ids)

    @staticmethod
    def parse_input(line, arrow=None, database=None, has_reaction_ids=False):
        reaction_id = None
        if has_reaction_ids:
            tokens = re.findall('(\w+)\s+(.*)', line.strip())[0]
            reaction_id = tokens[0]
            line = tokens[1]
        try:
            # This creates a Reaction instance for each reaction. each react is a sparse dict. also has cc
            reaction = Reaction.from_equation(line, arrow, database, reaction_id)
        except ParseException as e:
            logger.warning(str(e))
            reaction = Reaction({})
        return reaction

    @staticmethod
    def from_formulas(reaction_strings, database='bigg', arrow=None, has_reaction_ids=False, raise_exception=False):
        """
        parses a list of reactions in any format

        Arguments
        ---------
        reaction_strings : list
            A list of reactions in mol files.
        arrow : str
            The string used as the 'arrow' in each reaction (default: '<=>')

        Returns
        -------
        MetabolicModel
            A metabolic model
        """
        if arrow is None:
            arrow = '<=>'
        try:
            reactions = []
            not_balanced_count = 0
            for line in reaction_strings:

                reaction = MetabolicModel.parse_input(line, arrow, database, has_reaction_ids)

                if not reaction.is_balanced(fix_water=True, raise_exception=raise_exception):
                    not_balanced_count += 1
                    logger.warning('Model contains an unbalanced reaction: ' + line)
                    reaction = Reaction({})

                reactions.append(reaction)
                logger.debug('Adding reaction: ' + reaction.equation)

            if not_balanced_count > 0:
                warning_str = '%d out of the %d reactions are not chemically balanced' % (not_balanced_count,
                                                                                          len(reaction_strings))
                logger.debug(warning_str)
            return MetabolicModel.from_reactions(reactions, format, has_reaction_ids)

        except ValueError as e:
            if raise_exception:
                raise e
            else:
                logging.debug(str(e))
                return None

    @classmethod
    def from_kegg_formulas(cls, reaction_strings, arrow=None, has_reaction_ids=False, raise_exception=False):
        """
        Parses a list of reactions in KEGG format
        
        Arguments
        ---------
        reaction_strings : list
            A list of reactions in KEGG format.
        arrow : str
            The string used as the 'arrow' in each reaction (default: '<=>')
        has_reaction_ids : bool
            A flag indicating if there is a column of reaction IDs (separated from the reaction with whitespaces)
        raise_exception : bool
            If true, an exception will be raise if any reaction is unbalanced.

        Returns
        -------
        metabolic_model : MetabolicModel
        """
        if arrow is None:
            arrow = "<=>"

        try:
            reactions = []
            not_balanced_count = 0
            for line in reaction_strings:
                rid = None
                if has_reaction_ids:
                    tokens = re.findall('(\w+)\s+(.*)', line.strip())[0]
                    rid = tokens[0]
                    line = tokens[1]
                try:
                    # This creates a Reaction instance for each reaction. each react is a sparse dict. also has cc
                    reaction = Reaction.from_equation(line, arrow, rid)
                except ParseException as e:
                    logging.warning(str(e))
                    reaction = Reaction({})
                if not reaction.is_balanced(fix_water=True, raise_exception=raise_exception):
                    not_balanced_count += 1
                    logging.warning('Model contains an unbalanced reaction: ' + line)
                    reaction = Reaction({})
                reactions.append(reaction)
                logger.debug('Adding reaction: ' + reaction.equation)
            
            if not_balanced_count > 0:
                logger.debug('%d out of the %d reactions are not chemically balanced'
                             % (not_balanced_count, len(reaction_strings)))
            return cls.from_kegg_reactions(reactions, has_reaction_ids)
        
        except ValueError as e:
            if raise_exception:
                raise e
            else:
                logging.debug(str(e))
                return None

    def add_thermo(self, component_contribution_model):
        # check that all CIDs in the reaction are already cached by CC
        number_of_compounds, number_of_reactions = self.s_matrix.shape
        reactions = []
        for j in range(number_of_reactions):
            stoichiometry = {self.compound_ids[i]: self.s_matrix[i, j]
                             for i in range(number_of_compounds) if self.s_matrix[i, j] != 0}
            reaction = Reaction(stoichiometry)
            reactions.append(reaction)

        self.dg0, self.cov_dg0 = component_contribution_model.get_dg0_r_multi(reactions)

    def get_transformed_dg0(self, ph, ionic_strength, temperature):
        """
            returns the estimated dG0_prime and the standard deviation of
            each estimate (i.e. a measure for the uncertainty).
        """
        transformed_dg_component = self._get_transform_ddG0(ph=ph, ionic_strength=ionic_strength, temperature=temperature)
        dg0_prime = self.dg0 + transformed_dg_component
        dg0_std = np.matrix(np.sqrt(np.diag(self.cov_dg0))).T
        U, s, V = np.linalg.svd(self.cov_dg0, full_matrices=True)
        sqrt_sigma = np.matrix(U) * np.matrix(np.diag(s**0.5)) * np.matrix(V)
        return dg0_prime, dg0_std, sqrt_sigma

    def _get_transform_ddG0(self, ph, ionic_strength, temperature):
        """
        needed in order to calculate the transformed Gibbs energies of the 
        prediction_model reactions.
        
        Returns:
            an array (whose length is self.S.shape[1]) with the differences
            between DrG0_prime and DrG0. Therefore, one must add this array
            to the chemical Gibbs energies of reaction (DrG0) to get the 
            transformed values
        """
        cant_transform = []
        ddg0_compounds = np.matrix(np.zeros((self.s_matrix.shape[0], 1)))
        for i, cid in enumerate(self.compound_ids):
            compound = self.ccache.get_compound(cid)
            try:
                ddg0_compounds[i, 0] = compound.transform_pH7(ph, ionic_strength, temperature)
                if compound.inchi is None:
                    cant_transform.append(i)
            except KeyError:
                ddg0_compounds[i, 0] = 0
                cant_transform.append(i)
        
        ddg0_forward = np.dot(self.s_matrix.T, ddg0_compounds)

        # the dot product can be calculated even if some of the compounds couldn't be transformed.
        # here we find those reactions and nullify the transformation.
        indices = np.nonzero(self.s_matrix[cant_transform, :])
        skip_reactions = np.unique(indices[1])
        ddg0_forward[skip_reactions] = None

        return ddg0_forward
        
    def check_stoichiometric_balance(self, fix_water=False):
        elements, Ematrix = self.ccache.get_element_matrix(self.compound_ids)
        conserved = Ematrix.T * self.s_matrix

        if fix_water:
            # This part only looks for imbalanced oxygen and uses extra
            # H2O molecules (on either side of the reaction equation) to
            # balance them. Keep in mind that also the e- balance is affected
            # by the water (and hydrogen is not counted at all).
            if 'C00001' not in self.compound_ids:
                self.s_matrix = np.vstack([self.s_matrix, np.zeros((1, self.s_matrix.shape[1]))])
                self.compound_ids.append('C00001')
                elements, Ematrix = self.ccache.get_element_matrix(self.compound_ids)
            
            i_h2o = self.compound_ids.index('C00001')
            add_water = -conserved[elements.index('O'), :]
            self.s_matrix[i_h2o, :] += add_water
            conserved += Ematrix[i_h2o, :].T * add_water

        rxnFil = np.any(conserved[:, range(self.s_matrix.shape[1])], axis=0)
        unbalanced_ind = np.nonzero(rxnFil)[1]
        if len(unbalanced_ind) == 0:
            logging.warning('There are (%d) unbalanced reactions in S. Setting their coefficients to 0.' %
                            len(unbalanced_ind.flat))
            if self.rids is not None:
                logging.warning('These are the unbalanced reactions: ' +
                                ', '.join([self.rids[i] for i in unbalanced_ind.flat]))

            self.s_matrix[:, unbalanced_ind] = 0
        return self

    def write_reaction_by_index(self, r):
        stoichiometry = dict([(cid, self.s_matrix[i, r]) for i, cid in enumerate(self.compound_ids) if self.s_matrix[i, r] != 0])
        if self.rids is not None:
            reaction = Reaction(stoichiometry, reaction_id=self.rids[r])
        else:
            reaction = Reaction(stoichiometry)
        return reaction.equation
        
    def get_unidirectional_S(self):
        S_plus = np.copy(self.s_matrix)
        S_minus = np.copy(self.s_matrix)
        S_plus[self.s_matrix < 0] = 0
        S_minus[self.s_matrix > 0] = 0
        return S_minus, S_plus
