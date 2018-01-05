import csv
import logging

from component_contribution.prediction_model.model import ComponentContributionModel
from component_contribution.core.reaction import Reaction
from component_contribution.prediction_model.training_data import TrainingData
import numpy as np

logger = logging.getLogger('')
logger.setLevel(logging.INFO)


td = TrainingData()
cc = ComponentContributionModel.init()

G1 = np.matrix(cc.params['G1'])
G2 = np.matrix(cc.params['G2'])
G3 = np.matrix(cc.params['G3'])

############################################################################

reaction = Reaction.from_equation('C00002 + C00001 <=> C00008 + C00009')
print(reaction.equation)
# reaction = Reaction.from_equation('C00149 <=> C00122 + C00001')
file_name = 'reaction'

x, g = cc.decompose_reaction(reaction)
weights_rc = (x.T * G1).round(5)
weights_gc = (x.T * G2 + g.T * G3).round(5)
weights = weights_rc + weights_gc
orders = sorted(range(weights.shape[1]), key=lambda j: abs(weights[0, j]), reverse=True)

print(len(orders), len(td.dG0_prime), weights_gc.shape, weights_rc.shape)
print(weights.max(), weights.min())
output = csv.writer(open('%s_analysis.csv' % file_name, 'w'))
output.writerow(('Weight', 'dG\'0', 'dG0', 'reference', 'reaction'))

for i in orders:
    if abs(weights[0, i]) < 1e-7:
        continue

    a = weights_rc[0, i]
    a = td.dG0_prime[i]
    a = td.dG0[i]
    a = td.reference[i]
    a = td.description[i]
    output.writerow((weights_rc[0, i], td.dG0_prime[i], td.dG0[i], td.reference[i], td.description[i]))
