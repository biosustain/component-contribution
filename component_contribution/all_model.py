import re, csv, logging, copy
import numpy as np
from all_reaction import AllReaction
from kegg_reaction import KeggReaction
from compound_cacher import KeggCompoundCacher
from kegg_errors import KeggParseException
from CfB_functions import process_input_mets, _decompose_bigg_reaction, only_decompose

class AllModel(object):



    def __del__(self):
        self.ccache.dump()
    
    def __init__(self, S, input_metabolites, format, rids=None):
        self.S = S
        self.Smets = copy.copy(S)
        self.input_metabolites = input_metabolites
        self.rids = rids
        self.format = format

        assert len(self.input_metabolites) == self.S.shape[0]

        if self.rids is not None:
            assert len(self.rids) == self.S.shape[1]


        ccache = KeggCompoundCacher()

        # get compound data for all the input metabolites
        self.comp_data = process_input_mets(input_metabolites, ccache)

        self.cids = self.comp_data['output_ids']



        self.ccache = ccache

        # merge compound rows. this will leave a zero in transport reactions which this method currently can't handle!
        unique_cids = list(set(self.cids))
        for id in unique_cids:
            indeces = np.where(np.array(self.cids)==id)
            if indeces[0].size == 0:
                raise ValueError
            if indeces[0].size==1:
                continue
            else:

                self.S[indeces[0][0], :] += sum(self.S[indeces[0][1:], :])
                self.S = np.delete(self.S, (indeces[0][1:]), axis=0)
                self.cids = [i for j, i in enumerate(self.cids) if j not in indeces[0][1:]]

        # keep tab on transports
        transports=[]
        for p, i in enumerate(self.S.transpose()):
            if (i == 0).all():
                transports.append(p)
        self.transports = transports

        # DO THIS FOR OTHER FORMATS AS WELL
        # remove H+ from the stoichiometric matrix if it exists
        while 'C00080' in self.cids:
            i = self.cids.index('C00080')
            self.S = np.vstack((self.S[:i, :], self.S[i + 1:, :]))
            self.cids.pop(i)


    @staticmethod
    def from_file(fname, arrow='<=>', format='kegg', has_reaction_ids=False):
        """
        reads a file containing reactions in KEGG format
        
        Arguments:
           fname            - the filename to read
           arrow            - the string used as the 'arrow' in each reaction (default: '<=>')
           format           - the text file format provided ('kegg', 'tsv' or 'csv')
           has_reaction_ids - a boolean flag indicating if there is a column of
                              reaction IDs (separated from the reaction with
                              whitespaces)
        
        Return a KeggModel
        """
        fd = open(fname, 'r')
        if format == 'kegg':
            model = KeggModel.from_kegg_formulas(fd.readlines(), arrow, has_reaction_ids)
        elif format == 'tsv':
            model = KeggModel.from_csv(fd, has_reaction_ids=has_reaction_ids, delimiter='\t')
        elif format == 'csv':
            model = KeggModel.from_csv(fd, has_reaction_ids=has_reaction_ids, delimiter=None)
        fd.close()
        return model
    
    @staticmethod
    def from_csv(fd, has_reaction_ids=True, delimiter=None):
        csv_reader = csv.reader(fd, delimiter=delimiter)
        if has_reaction_ids:
            rids = csv_reader.next()
            rids = rids[1:]
        else:
            rids = None
        S = []
        cids = []
        for i, row in enumerate(csv_reader):
            cids.append(row[0])
            S.append([float(x) for x in row[1:]])
        S = np.array(S)

        return KeggModel(S, cids, rids)

    @staticmethod
    def from_reactions(reactions, format, has_reaction_ids=False):
        if has_reaction_ids:
            rids = [r.rid for r in reactions]
        else:
            rids = None

        input_ids = set()
        for reaction in reactions:
            input_ids = input_ids.union(reaction.keys())

        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        input_ids = sorted(input_ids)
        S = np.matrix(np.zeros((len(input_ids), len(reactions))))
        for i, reaction in enumerate(reactions):
            S[:, i] = np.matrix(reaction.dense(input_ids))

        logging.debug('Successfully loaded %d reactions (involving %d unique compounds)' %
                      (S.shape[1], S.shape[0]))
        return AllModel(S, input_ids, format, rids)

    @staticmethod
    def from_kegg_reactions(kegg_reactions, has_reaction_ids=False):
        if has_reaction_ids:
            rids = [r.rid for r in kegg_reactions]
        else:
            rids = None

        cids = set()
        for reaction in kegg_reactions:
            cids = cids.union(reaction.keys())
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        cids = sorted(cids)
        S = np.matrix(np.zeros((len(cids), len(kegg_reactions))))
        for i, reaction in enumerate(kegg_reactions):
            S[:, i] = np.matrix(reaction.dense(cids))
        
        logging.debug('Successfully loaded %d reactions (involving %d unique compounds)' %
                      (S.shape[1], S.shape[0]))
        return KeggModel(S, cids, rids)

    @staticmethod
    def parse_input(line, format, arrow='<=>', has_reaction_ids = False):
        rid = None
        if has_reaction_ids:
            tokens = re.findall('(\w+)\s+(.*)', line.strip())[0]
            rid = tokens[0]
            line = tokens[1]
        try:
            # This creates a KeggReaction instance for each reaction. each react is a sparse dict. also has cc
            reaction = AllReaction.parse_formula(line, format, arrow, rid)
        except KeggParseException as e:
            logging.warning(str(e))
            reaction = AllReaction({})
        return reaction

    @staticmethod
    def from_formulas(reaction_strings, format='bigg', arrow='<=>', has_reaction_ids=False,
                      raise_exception=False):
        """
        parses a list of reactions in any format

        Arguments:
           reaction_strings - a list of reactions in mol files
           arrow            - the string used as the 'arrow' in each reaction (default: '<=>')


        Return values:
           S     - a stoichiometric matrix
           cids  - the KEGG compound IDs in the same order as the rows of S
        """
        try:
            reactions = []
            not_balanced_count = 0
            for line in reaction_strings:

                reaction = AllModel.parse_input(line, format, arrow, has_reaction_ids)

                # uncomment balance check!
                          ############
                if False: # not reaction.is_balanced(fix_water=True, raise_exception=raise_exception):
                          ############
                    not_balanced_count += 1
                    logging.warning('Model contains an unbalanced reaction: ' + line)
                    reaction = KeggReaction({})
                reactions.append(reaction)
                logging.debug('Adding reaction: ' + reaction.write_formula())

            if not_balanced_count > 0:
                warning_str = '%d out of the %d reactions are not chemically balanced' % \
                              (not_balanced_count, len(reaction_strings))
                logging.debug(warning_str)
            return AllModel.from_reactions(reactions, format, has_reaction_ids)

        except ValueError as e:
            if raise_exception:
                raise e
            else:
                logging.debug(str(e))
                return None

    @staticmethod
    def from_kegg_formulas(reaction_strings, arrow='<=>', has_reaction_ids=False,
                           raise_exception=False):
        """
        parses a list of reactions in KEGG format
        
        Arguments:
           reaction_strings - a list of reactions in KEGG format
           arrow            - the string used as the 'arrow' in each reaction (default: '<=>')
           has_reaction_ids - a boolean flag indicating if there is a column of
                              reaction IDs (separated from the reaction with
                              whitespaces)
        
        Return values:
           S     - a stoichiometric matrix
           cids  - the KEGG compound IDs in the same order as the rows of S
        """
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
                    # This creates a KeggReaction instance for each reaction. each react is a sparse dict. also has cc
                    reaction = KeggReaction.parse_formula(line, arrow, rid)
                except KeggParseException as e:
                    logging.warning(str(e))
                    reaction = KeggReaction({})
                if not reaction.is_balanced(fix_water=True, raise_exception=raise_exception):
                    not_balanced_count += 1
                    logging.warning('Model contains an unbalanced reaction: ' + line)
                    reaction = KeggReaction({})
                reactions.append(reaction)
                logging.debug('Adding reaction: ' + reaction.write_formula())
            
            if not_balanced_count > 0:
                warning_str = '%d out of the %d reactions are not chemically balanced' % \
                              (not_balanced_count, len(reaction_strings))
                logging.debug(warning_str)
            return KeggModel.from_kegg_reactions(reactions, has_reaction_ids)
        
        except ValueError as e:
            if raise_exception:
                raise e
            else:
                logging.debug(str(e))
                return None

    def add_thermo(self, cc):

        # decompose the separate compnents from the different ids

        # doing for kegg

        S_kegg = self.comp_data['S_kegg']

        bigg_to_kegg = self.comp_data['bigg_to_kegg']
        bigg_to_kegg_dict = self.comp_data['bigg_to_kegg_dict']
        kegg_ids = [bigg_to_kegg_dict[id] for id in bigg_to_kegg]

        Nc, Nr = S_kegg.shape
        kegg_reactions = []
        for j in xrange(Nr):
            sparse = {kegg_ids[i]:S_kegg[i,j] for i in xrange(Nc) if S_kegg[i,j] != 0}
            reaction = KeggReaction(sparse)
            kegg_reactions.append(reaction)

        X_kegg, G_kegg = only_decompose(cc, kegg_reactions)



        self.dG0, self.cov_dG0 = cc.get_dG0_r_multi(reactions)




        # do for smile
        try:
            if True:
                x, g, model_ids_to_replace_per_reac = _decompose_bigg_reaction(self, reaction, bigg_dict)

            else:
                x, g = self._decompose_reaction(reaction)

        except inchi2gv.GroupDecompositionError:
            pass



        # do for bigg non kegg



        self.dG0, self.cov_dG0, self.model_ids_to_replace = cc.get_dG0_r_multi(reactions)
        
    def get_transformed_dG0(self, pH, I, T):
        """
            returns the estimated dG0_prime and the standard deviation of
            each estimate (i.e. a measure for the uncertainty).
        """
        self.transformed_dG_component = self._get_transform_ddG0(pH=pH, I=I, T=T)
        self.dG0_prime = self.dG0 + self.transformed_dG_component
        dG0_std = np.matrix(np.sqrt(np.diag(self.cov_dG0))).T
        U, s, V = np.linalg.svd(self.cov_dG0, full_matrices=True)
        sqrt_Sigma = np.matrix(U) * np.matrix(np.diag(s**0.5)) * np.matrix(V)
        return self.dG0_prime, dG0_std, sqrt_Sigma

    def _get_transform_ddG0(self, pH, I, T):
        """
        needed in order to calculate the transformed Gibbs energies of the 
        model reactions.
        
        Returns:
            an array (whose length is self.S.shape[1]) with the differences
            between DrG0_prime and DrG0. Therefore, one must add this array
            to the chemical Gibbs energies of reaction (DrG0) to get the 
            transformed values
        """
        cant_transform = []
        ddG0_compounds = np.matrix(np.zeros((self.S.shape[0], 1)))
        for i, cid in enumerate(self.cids):
            if cid in self.comp_data['all_comp_data']:
                comp = self.comp_data['all_comp_data'][cid]
            else:
                comp = None
            #    comp = self.ccache.get_compound(cid)

            try:
                ddG0_compounds[i, 0] = comp.transform_pH7(pH, I, T)
                if comp.inchi == None: cant_transform.append(i)
            except:
                ddG0_compounds[i, 0] = 0
                cant_transform.append(i)
        
        ddG0_forward = np.dot(self.S.T, ddG0_compounds)

        # the dot product can be calculated even if some of the compounds couldn't be transformed.
        # here we find those reactions and nullify the transformation.
        inds = np.nonzero(self.S[cant_transform, :])
        reacs_that_shouldnt_be_transformed = np.unique(inds[1])
        ddG0_forward[reacs_that_shouldnt_be_transformed] = None

        return ddG0_forward
        
    def check_S_balance(self, fix_water=False):
        elements, Ematrix = self.ccache.get_element_matrix(self.cids)
        conserved = Ematrix.T * self.S

        if fix_water:
            # This part only looks for imbalanced oxygen and uses extra
            # H2O molecules (on either side of the reaction equation) to
            # balance them. Keep in mind that also the e- balance is affected
            # by the water (and hydrogen is not counted at all).
            if 'C00001' not in self.cids:
                self.S = np.vstack([self.S, np.zeros((1, self.S.shape[1]))])
                self.cids.append('C00001')
                elements, Ematrix = self.ccache.get_element_matrix(self.cids)
            
            i_h2o = self.cids.index('C00001')
            add_water = -conserved[elements.index('O'), :]
            self.S[i_h2o, :] += add_water
            conserved += Ematrix[i_h2o, :].T * add_water

        rxnFil = np.any(conserved[:,range(self.S.shape[1])],axis=0)
        unbalanced_ind = np.nonzero(rxnFil)[1]
        if unbalanced_ind != []:
            logging.warning('There are (%d) unbalanced reactions in S. ' 
                            'Setting their coefficients to 0.' % 
                            len(unbalanced_ind.flat))
            if self.rids is not None:
                logging.warning('These are the unbalanced reactions: ' +
                                ', '.join([self.rids[i] for i in unbalanced_ind.flat]))
                    
            self.S[:, unbalanced_ind] = 0
        return self

    def write_reaction_by_index(self, r):
        sparse = dict([(cid, self.S[i, r]) for i, cid in enumerate(self.cids)
                       if self.S[i, r] != 0])
        if self.rids is not None:
            reaction = KeggReaction(sparse, rid=self.rids[r])
        else:
            reaction = KeggReaction(sparse)
        return reaction.write_formula()
        
    def get_unidirectional_S(self):
        S_plus = np.copy(self.S)
        S_minus = np.copy(self.S)
        S_plus[self.S < 0] = 0
        S_minus[self.S > 0] = 0
        return S_minus, S_plus
        