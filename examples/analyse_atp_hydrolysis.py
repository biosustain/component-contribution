import logging, csv
import numpy as np
logger = logging.getLogger('')
logger.setLevel(logging.INFO)
from component_contribution.component_contribution_trainer import ComponentContribution
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.training_data import TrainingData

def main():

    td = TrainingData()
    cc = ComponentContribution.init()

    G1 = np.matrix(cc.params['G1'])
    G2 = np.matrix(cc.params['G2'])
    G3 = np.matrix(cc.params['G3'])

    ############################################################################

    #reaction = KeggReaction.parse_formula('C00002 + C00001 <=> C00008 + C00009'); fname = 'atpase';
    reaction = KeggReaction.parse_formula('C00149 <=> C00122 + C00001'); fname = 'empty';

    x, g = cc._decompose_reaction(reaction)
    weights_rc = (x.T * G1).round(5)
    weights_gc = (x.T * G2 + g.T * G3).round(5)
    weights = weights_rc + weights_gc
    orders = sorted(range(weights.shape[1]),
                    key=lambda j:abs(weights[0, j]), reverse=True)

    output = csv.writer(open('%s_analysis.csv' % fname, 'w'))
    output.writerow(('Weight', 'dG\'0', 'dG0', 'reference', 'reaction'))
    for j in orders:
        if abs(weights[0, j]) < 1e-7:
            continue
                              
        output.writerow((weights_rc[0, j], td.dG0_prime[j], td.dG0[j],
                         td.reference[j], td.description[j]))

if __name__ == '__main__':
    main()