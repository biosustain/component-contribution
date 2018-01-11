import logging
import re

import numpy as np

from component_contribution.compound_cache import CompoundCache
from component_contribution.exceptions import ParseException

logger = logging.getLogger(__name__)


class Reaction(object):
    """
    Reaction representation.

    Attributes
    ----------
    stoichiometry : dict
        The stoichiometric representation of the reaction.
    arrow : str
        The arrow token to split the reaction.
    reaction_id : str
        The identifier of the reaction.

    """
    def __init__(self, stoichiometry, arrow=None, database=None, reaction_id=None):
        if arrow is None:
            arrow = '<=>'

        for compound_id, coefficient in stoichiometry.items():
            if not isinstance(coefficient, (float, int)):
                raise ValueError('All stoichiometric coefficients in Reaction must be integers or floats')

        self.compound_cache = CompoundCache.get_instance()

        # translate compound ids
        self.stoichiometry = {}
        for compound_id, coefficient in stoichiometry.items():
            if coefficient != 0:
                compound = self.compound_cache.get_compound(compound_id)
                self.stoichiometry[compound.id] = coefficient

        self.arrow = arrow
        self.database = database
        self.reaction_id = reaction_id


    @staticmethod
    def __format_compound(compound_id, coefficient):
        if coefficient == 1:
            return compound_id
        else:
            return "%s %s" % (coefficient, compound_id)

    @classmethod
    def from_equation(cls, equation, arrow='<=>', database=None, reaction_id=None):
        """
        Parses a two-sided formula such as: 2 C00001 => C00002 + C00003

        Returns
        -------

            The set of substrates, products and the direction of the reaction
        """
        tokens = equation.split(arrow)
        if len(tokens) < 2:
            raise ParseException('Reaction does not contain the arrow sign (%s): %s' % (arrow, equation))
        if len(tokens) > 2:
            raise ParseException('Reaction contains more than one arrow sign (%s): %s' % (arrow, equation))

        left = tokens[0].strip()
        right = tokens[1].strip()

        sparse_reaction = {}
        for met_id, count in cls.parse_reaction_formula_side(left).items():
            sparse_reaction[met_id] = sparse_reaction.get(met_id, 0) - count

        for met_id, count in cls.parse_reaction_formula_side(right).items():
            sparse_reaction[met_id] = sparse_reaction.get(met_id, 0) + count

        return cls(sparse_reaction, arrow, database=database, reaction_id=reaction_id)

    @staticmethod
    def parse_reaction_formula_side(s):
        """
        Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        Ignores stoichiometry.

        Returns:
            The set of CIDs.
        """
        if s.strip() == "null":
            return {}

        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)
            if len(tokens) == 0:
                continue
            if len(tokens) == 1:
                amount = 1
                key = member
            else:
                try:
                    amount = float(tokens[0])
                except ValueError:
                    raise ParseException("Non-specific reaction: %s" % s)
                key = tokens[1]

            try:
                compound_bag[key] = compound_bag.get(key, 0) + amount
            except ValueError:
                raise ParseException("Non-specific reaction: %s" % s)

        return compound_bag

    @property
    def compound_ids(self):
        return set(self.stoichiometry.keys())

    @property
    def quotient(self):
        products = sum(self.compound_cache.get_compound(compound_id).concentration ** coefficient
                       for compound_id, coefficient in self.stoichiometry.items() if coefficient > 0)

        substrates = sum(self.compound_cache.get_compound(compound_id).concentration ** coefficient
                         for compound_id, coefficient in self.stoichiometry.items() if coefficient < 0)

        return products/substrates

    @property
    def equation(self):
        """
        Returns a string representing the chemical conversion.

        Returns
        -------
        str
        """
        left = []
        right = []
        for compound_id, coefficient in sorted(self.stoichiometry.items()):
            if coefficient < 0:
                left.append(self.__format_compound(compound_id, -coefficient))
            elif coefficient > 0:
                right.append(self.__format_compound(compound_id, coefficient))
        return "%s %s %s" % (' + '.join(left), self.arrow, ' + '.join(right))

    def transformed_gibbs_energy(self, ph, ionic_strength, temperature):
        """
        Needed in order to calculate the transformed Gibbs energies of reactions.

        The difference between DrG0_prime and DrG0 for this reaction.
        Therefore, this value must be added to the chemical Gibbs energy
        of reaction (DrG0) to get the transformed value.

        Parameters
        ----------
        ph : float
            The pH value [0..14].
        ionic_strength : float
            The ionic strength.
        temperature : float
            The temperature in Kelvin.

        Returns
        -------
        float
            The transformed energy
        """
        if not 0 <= ph <= 14:
            raise ValueError("ph must be between 0 and 14")

        change_forward_gibbs_energy = 0
        for compound_id, coefficient in self.stoichiometry.items():
            compound = self.compound_cache.get_compound(compound_id)
            change_forward_gibbs_energy += coefficient * compound.transform_ph7(ph, ionic_strength, temperature)
        return change_forward_gibbs_energy

    def atom_bag(self, raise_exception=False):
        """
        Use for checking for mass balance.

        An atom_bag of the differences between the sides of the reaction. E.g. if there is one extra C on the
        left-hand side, the result will be {'C': -1}.

        Returns
        -------
        atom_bag : dict
            The count of elements balance.
        """
        try:
            compound_ids = list(self.stoichiometry.keys())
            coefficients = map(self.stoichiometry.__getitem__, compound_ids)
            coefficients = np.array(list(coefficients))

            cached_compound_ids = set(self.compound_cache.compound_ids)
            if not cached_compound_ids.issuperset(compound_ids):
                missing_compound_ids = set(compound_ids).difference(cached_compound_ids)
                error_message = 'The following compound IDs are not in the cache,  make sure they appear ' \
                                'in kegg_additions.tsv and then run compound_cache.py: ' + \
                                ', '.join(sorted(missing_compound_ids))
                raise ValueError(error_message)

            elements, element_matrix = self.compound_cache.get_element_matrix(compound_ids)
            conserved = coefficients * element_matrix

            if np.any(np.isnan(conserved), 1):
                error_message = 'cannot test reaction balancing because of unspecific ' + \
                              'compound formulas: %s' % self.equation
                raise ValueError(error_message)

            atom_bag = {}
            if np.any(conserved != 0, axis=1):
                logger.debug('unbalanced reaction: %s' % self.equation)
                for j, c in enumerate(conserved.flat):
                    if c != 0:
                        logger.debug('there are %d more %s atoms on the right-hand side' % (c, elements[j]))
                        atom_bag[str(elements[j])] = c
            return atom_bag

        except ValueError:
            if raise_exception:
                raise ValueError
            else:
                logging.warning(str(ValueError('could not balance reaction' + str(self))))
                return None

    def is_balanced(self, fix_water=False, raise_exception=False):
        reaction_atom_bag = self.atom_bag(raise_exception)

        if reaction_atom_bag is None:  # this means some compound formulas are missing
            return False

        if fix_water and 'O' in reaction_atom_bag:
            self.stoichiometry.setdefault('C00001', 0)
            self.stoichiometry['C00001'] += -reaction_atom_bag['O']
            if self.stoichiometry['C00001'] == 0:
                del self.stoichiometry['C00001']
            reaction_atom_bag = self.atom_bag()

        return len(reaction_atom_bag) == 0

    def is_empty(self):
        return len(self.stoichiometry) == 0

    def dense(self, ids):
        s = np.matrix(np.zeros((len(ids), 1)))
        for compound_id, coefficient in self.stoichiometry.items():
            s[ids.index(compound_id), 0] = coefficient
        return s

    def transform_ddg0(self, ph, ionic_strength, temperature):
        """
        Needed in order to calculate the transformed Gibbs energies of
        reactions.

        Returns:
            The difference between DrG0_prime and DrG0 for this reaction.
            Therefore, this value must be added to the chemical Gibbs
            energy of reaction (DrG0) to get the transformed value.
        """
        ddg0_forward = 0
        for compound_id, coefficient in self.stoichiometry.items():
            comp = self.compound_cache.get_compound(compound_id)
            ddg0_forward += coefficient * comp.transform_ph7(ph, ionic_strength, temperature)
        return ddg0_forward

