import os
import logging
import csv
import numpy as np
from scipy.io import savemat
from component_contribution.thermodynamics.constants import R, F
from component_contribution.compound_cache import CompoundCache
from component_contribution.core.reaction import Reaction

__DIR = os.path.realpath(os.path.dirname(__file__))
__DATA_DIR = os.path.join(__DIR, "..", "..", 'data')

TECRDB_PATH = os.path.join(__DATA_DIR, "TECRDB.tsv")
TECRDB_WEIGHT = 1.0

FORMATION_PATH = os.path.join(__DATA_DIR, "formation_energies_transformed.tsv")
FORMATION_WEIGHT = 1.0

REDOX_PATH = os.path.join(__DATA_DIR, "redox.tsv")
REDOX_WEIGHT = 1.0


class TrainingData(object):

    def __init__(self):
        self.ccache = CompoundCache.get_instance()

        # this loads all the training data and returns is as a dict in thermo_params.
        # it also defines which metabolites not to decompose, there are some metabolites for which there are formation
        # energies, but they are still decomposed.
        thermo_params, self.cids_that_dont_decompose = TrainingData.get_all_thermo_params()

        # cids is a set of all the cids used in the training data this is saved as cids_joined outside
        # this function (inside ComponentContribution() )
        cids = set()
        for d in thermo_params:
            cids = cids.union(d['reaction'].compound_ids)
        cids = sorted(cids)
        
        # convert the list of reactions in sparse notation into a full
        # stoichiometric matrix, where the rows (compounds) are according to the
        # CID list 'cids'.
        self.S = np.zeros((len(cids), len(thermo_params)))
        for k, a in enumerate(thermo_params):
            for compound, coefficient in a['reaction'].metabolites.items():
                self.S[cids.index(compound.id), k] = coefficient
            
        self.cids = cids

        self.dG0_prime = np.array([d['dG\'0'] for d in thermo_params])
        self.T = np.array([d['T'] for d in thermo_params])
        self.I = np.array([d['I'] for d in thermo_params])
        self.pH = np.array([d['pH'] for d in thermo_params])
        self.pMg = np.array([d['pMg'] for d in thermo_params])
        self.weight = np.array([d['weight'] for d in thermo_params])
        self.reference = [d['reference'] for d in thermo_params]
        self.description = [d['description'] for d in thermo_params]
        rxn_inds_to_balance = [i for i in range(len(thermo_params)) if thermo_params[i]['balance']]

        # train only on those reactions that are elementally balanced. if O is missing, adds water.
        self.balance_reactions(rxn_inds_to_balance)

        # this corrects for the non - standard condition in which the experimental data where measured.
        # It back calculates dG0 saved in self.dG0
        self.dG0 = None
        self.reverse_transform()

    def __del__(self):
        self.ccache.dump()

    def to_mat(self, file_name):
        """
        Write all training data to a Matlab file.
            
        Arguments
        ---------
        file_name : str
            Where the data will be written.
        """
        d = {'dG0_prime': self.dG0_prime,
             'dG0': self.dG0,
             'T': self.T,
             'I': self.I,
             'pH': self.pH,
             'pMg': self.pMg,
             'weight': self.weight,
             'cids': self.cids,
             'S': self.S}
        savemat(file_name, d, oned_as='row')

    def to_csv(self, fname):
        csv_output = csv.writer(open(fname, 'w'))
        csv_output.writerow(['reaction', 'T', 'I', 'pH', 'reference', 'dG0', 'dG0_prime'])
        for j in range(self.S.shape[1]):
            stoichiometry = {self.ccache.get_compound(self.cids[i]): self.S[i, j] for i in range(self.S.shape[0])}
            reaction = Reaction()
            reaction.add_metabolites(stoichiometry)
            csv_output.writerow([reaction.reaction, self.T[j], self.I[j], self.pH[j],
                                 self.reference[j], self.dG0[j], self.dG0_prime[j]])

    @staticmethod
    def str2double(s):
        """
            casts a string to float, but if the string is empty return NaN
        """
        if s == '':
            return np.nan
        else:
            return float(s)

    @staticmethod
    def read_tecrdb(fname, weight):
        """Read the raw data of TECRDB (NIST)"""
        thermo_params = []  # columns are: reaction, dG'0, T, I, pH, pMg, weight, balance?

        headers = ["URL", "REF_ID", "METHOD", "EVAL", "EC", "ENZYME NAME", "REACTION IN KEGG IDS",
                   "REACTION IN COMPOUND NAMES", "K", "K'", "T", "I", "pH", "pMg"]

        for row_list in csv.reader(open(fname, 'r'), delimiter='\t'):
            if isinstance(row_list, list) and len(row_list) == 0:
                continue
            row = dict(zip(headers, row_list))
            if (row['K\''] == '') or (row['T'] == '') or (row['pH'] == ''):
                continue
            
            # parse the reaction
            reaction = Reaction.from_equation(row['REACTION IN KEGG IDS'], arrow='=')

            # calculate dG'0
            dG0_prime = -R * TrainingData.str2double(row['T']) * np.log(TrainingData.str2double(row['K\'']))
            try:
                thermo_params.append({'reaction': reaction,
                                      'dG\'0' : dG0_prime,
                                      'T': TrainingData.str2double(row['T']), 
                                      'I': TrainingData.str2double(row['I']),
                                      'pH': TrainingData.str2double(row['pH']),
                                      'pMg': TrainingData.str2double(row['pMg']),
                                      'weight': weight,
                                      'balance': True,
                                      'reference': row['REF_ID'],
                                      'description': row['REACTION IN COMPOUND NAMES']})
            except ValueError:
                raise Exception('Cannot parse row: ' + str(row))

        logging.debug('Successfully added %d reactions from TECRDB' % len(thermo_params))
        return thermo_params
        
    @staticmethod
    def read_formations(fname, weight):
        """Read the Formation Energy data"""
        
        # columns are: reaction, dG'0, T, I, pH, pMg, weight, balance?
        thermo_params = []
        cids_that_dont_decompose = set()
        
        # fields are: cid, name, dG'0, pH, I, pMg, T, decompose?,
        #             compound_ref, remark
        for row in csv.DictReader(open(fname, 'r'), delimiter='\t'):
            if int(row['decompose']) == 0:
                cids_that_dont_decompose.add(row['cid'])
            if row['dG\'0'] != '':
                rxn = Reaction({row['cid']: 1})
                thermo_params.append({'reaction': rxn,
                                      'dG\'0': TrainingData.str2double(row['dG\'0']),
                                      'T': TrainingData.str2double(row['T']), 
                                      'I': TrainingData.str2double(row['I']),
                                      'pH': TrainingData.str2double(row['pH']),
                                      'pMg': TrainingData.str2double(row['pMg']),
                                      'weight': weight,
                                      'balance': False,
                                      'reference': row['compound_ref'],
                                      'description': row['name'] + ' formation'})

        logging.debug('Successfully added %d formation energies' % len(thermo_params))
        return thermo_params, cids_that_dont_decompose
        
    @staticmethod
    def read_redox(fname, weight):
        """Read the Reduction potential data"""
        # columns are: reaction, dG'0, T, I, pH, pMg, weight, balance?
        thermo_params = []

        # fields are: name, CID_ox, nH_ox, charge_ox, CID_red,
        #             nH_red, charge_red, E'0, pH, I, pMg, T, ref
        for row in csv.DictReader(open(fname, 'r'), delimiter='\t'):
            delta_nH = TrainingData.str2double(row['nH_red']) - \
                       TrainingData.str2double(row['nH_ox'])
            delta_charge = TrainingData.str2double(row['charge_red']) - TrainingData.str2double(row['charge_ox'])
            delta_e = delta_nH - delta_charge
            dG0_prime = -F * TrainingData.str2double(row['E\'0']) * delta_e
            rxn = Reaction({row['CID_ox']: -1, row['CID_red']: 1})
            thermo_params.append({'reaction':    rxn,
                                  'dG\'0':       dG0_prime,
                                  'T':           TrainingData.str2double(row['T']),
                                  'I':           TrainingData.str2double(row['I']),
                                  'pH':          TrainingData.str2double(row['pH']),
                                  'pMg':         TrainingData.str2double(row['pMg']),
                                  'weight':      weight,
                                  'balance':     False,
                                  'reference':   row['ref'],
                                  'description': row['name'] + ' redox'})

        logging.debug('Successfully added %d redox potentials' % len(thermo_params))
        return thermo_params
    
    @staticmethod
    def get_all_thermo_params():
        # this loads all the training data and returns is as a dict in thermo_params.
        # it also defines which metabolites not to decompose, there are some metabolites for which there are formation
        # energies, but they are still decomposed.
        tecrdb_params = TrainingData.read_tecrdb(TECRDB_PATH, TECRDB_WEIGHT)
        formation_params, cids_that_dont_decompose = TrainingData.read_formations(FORMATION_PATH, FORMATION_WEIGHT)
        redox_params = TrainingData.read_redox(REDOX_PATH, REDOX_WEIGHT)
        
        thermo_params = tecrdb_params + formation_params + redox_params
        return thermo_params, cids_that_dont_decompose

    @property
    def water_id(self):
        return self.ccache.compound_id2inchi_key["C00001"]

    @property
    def proton_id(self):
        return self.ccache.compound_id2inchi_key['C00080']
    
    def balance_reactions(self, reaction_indices_to_balance):
        """
        use the chemical formulas from the InChIs to verify that each and every
        reaction is balanced
        """
        elements, element_matrix = self.ccache.get_element_matrix(self.cids)
        compound_indices_without_formula = list(np.nonzero(np.any(np.isnan(element_matrix), 1))[0].flat)
        element_matrix[np.isnan(element_matrix)] = 0

        S_without_formula = self.S[compound_indices_without_formula, :]
        reaction_indices_without_formula = np.nonzero(np.any(S_without_formula != 0, 0))[0]
        reaction_indices_to_balance = set(reaction_indices_to_balance).difference(reaction_indices_without_formula)

        # need to check that all elements are balanced (except H, but including e-)
        # if only O is not balanced, add water molecules
        if 'O' in elements:
            i_h2o = self.cids.index(self.water_id)
            j_o = elements.index('O')
            conserved = np.dot(element_matrix.T, self.S)
            for k in reaction_indices_to_balance:
                self.S[i_h2o, k] = self.S[i_h2o, k] - conserved[j_o, k]

        # recalculate conservation matrix
        conserved = element_matrix.T * self.S
        
        reaction_indices_to_remove = [k for k in reaction_indices_to_balance if np.any(conserved[:, k] != 0, 0)]

        # create an output in the logging file. is all this indexing necessary?
        for k in reaction_indices_to_remove:
            stoichiometry = {}
            for i in np.nonzero(self.S[:, k])[0]:
                stoichiometry[self.ccache.get_compound(self.cids[i])] = self.S[i, k]
            reaction = Reaction()
            reaction.add_metabolites(stoichiometry)
            logging.debug('unbalanced reaction #%d: %s' % (k, reaction.reaction))
            for j in np.where(conserved[:, k])[0].flat:
                logging.debug('there are %d more %s atoms on the right-hand side' % (conserved[j, k], elements[j]))
        
        reaction_indices_to_keep = set(range(self.S.shape[1])).difference(reaction_indices_to_remove)
        
        reaction_indices_to_keep = sorted(reaction_indices_to_keep)
        
        self.S = self.S[:, reaction_indices_to_keep]
        self.dG0_prime = self.dG0_prime[reaction_indices_to_keep]
        self.T = self.T[reaction_indices_to_keep]
        self.I = self.I[reaction_indices_to_keep]
        self.pH = self.pH[reaction_indices_to_keep]
        self.pMg = self.pMg[reaction_indices_to_keep]
        self.weight = self.weight[reaction_indices_to_keep]
        self.reference = [self.reference[i] for i in reaction_indices_to_keep]
        self.description = [self.description[i] for i in reaction_indices_to_keep]

        logging.debug('After removing %d unbalanced reactions, the stoichiometric matrix contains: '
                      '%d compounds and %d reactions' % (len(reaction_indices_to_remove), self.S.shape[0], self.S.shape[1]))

    def reverse_transform(self):
        """
        Calculate the reverse transform for all reactions in training_data.
        the training data gives us dG0', to reverse transform each reaction (i) we correct each individual
        compound (j) and then we add them
        """
        number_of_reactions = self.S.shape[1]
        reverse_ddG0 = np.zeros(number_of_reactions)
        self.I[np.isnan(self.I)] = 0.25  # default ionic strength is 0.25M
        self.pMg[np.isnan(self.pMg)] = 14  # default pMg is 14
        for i in range(number_of_reactions):
            for j in np.nonzero(self.S[:, i])[0]:
                cid = self.cids[j]
                if cid == self.proton_id:  # H+ should be ignored in the Legendre transform
                    continue
                comp = self.ccache.get_compound(cid)
                ddG0 = comp.transform_ph7(self.pH[i], self.I[i], self.T[i])
                reverse_ddG0[i] = reverse_ddG0[i] + ddG0 * self.S[j, i]

        self.dG0 = self.dG0_prime - reverse_ddG0


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=
        'Prepare all thermodynamic training data in a .mat file for running '
        'component contribution.')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                       help='the path to the .mat file that should be written '
                       'containing the training data')
    
    args = parser.parse_args()
    td = TrainingData()
    td.to_mat(args.outfile)
