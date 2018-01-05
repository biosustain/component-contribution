import csv
import os

import numpy as np

from component_contribution.core.metabolic_model import MetabolicModel
# import matplotlib.pyplot as plt
from component_contribution.prediction_model.model import ComponentContributionModel

REPORT_CACHE_FNAME = '../cache/report_gc.csv'
cid2dG0 = {}
if not os.path.exists(REPORT_CACHE_FNAME):
    fp = open(REPORT_CACHE_FNAME, 'w')
    cc = ComponentContributionModel()
    cc.train()
    csv_out = csv.writer(fp)
    csv_out.writerow(['cid', 'dG0_f'])
    for compound_id in cc.ccache.compound_ids:
        dG0_f = cc.get_major_ms_dG0_f(compound_id)
        csv_out.writerow([compound_id, '%8.2f' % dG0_f])
        cid2dG0[compound_id] = dG0_f
else:
    for row in csv.DictReader(open(REPORT_CACHE_FNAME, 'r')):
        cid2dG0[row['cid']] = float(row['dG0_f'])

REACTION_FNAME = 'tests/report_gc_reactions.txt'
reaction_strings = open(REACTION_FNAME, 'r').readlines()
model = MetabolicModel.from_kegg_formulas(reaction_strings)

# compare the dG0_r of the prediction_model to the one we get if multiplying
# the prediction_model stoichiometric matrix by the formation energies
model_dG0_f = np.matrix([cid2dG0[cid] for cid in model.compound_ids]).T
model_dG0_r = model.s_matrix.T * model_dG0_f

cc2 = ComponentContributionModel()
model.add_thermo(cc2)
plt.plot(model.dG0, model_dG0_r, '.')
plt.show()

#dG0_prime, dG0_std = prediction_model.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)

