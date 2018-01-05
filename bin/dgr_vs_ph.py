
import sys

import numpy as np
from scipy.io import savemat

from component_contribution.io.kegg_parser import KEGGData
from component_contribution.core.metabolic_model import MetabolicModel
from component_contribution.prediction_model.model import ComponentContributionModel


def kegg_file_to_models(pathway_file):
    kegg_file = KEGGData.from_kegg_file(pathway_file)
    entries = kegg_file.entries()
    pathways = []
    for entry in entries:
        fields = kegg_file[entry]
        rids, fluxes, reactions = KEGGData.parse_reaction_module(fields)
        bounds = KEGGData.parse_bound_module(fields)
        model = MetabolicModel.from_kegg_formulas(reactions)
        model.rids = rids
        pH = fields.get_float_field('PH', 7.5)
        I = fields.get_float_field('I', 0.2)
        T = fields.get_float_field('T', 298.15)
        pathways.append({'entry': entry, 'prediction_model': model, 'fluxes': fluxes,
                         'bounds': bounds, 'pH': pH, 'I': I, 'T': T})
    return pathways


if __name__ == '__main__':
    file_name = sys.argv[1]
    
    pathways = kegg_file_to_models(file_name)

    pH_range = np.arange(5, 9.0001, 0.1)
    I = 0.2
    T = 298.15

    cc = ComponentContributionModel.init()

    cc_dict = {'pH': np.matrix(pH_range).T}
    for p in pathways:
        entry = p['entry']
        p['prediction_model'].add_thermo(cc)
        
        dG0_primes = []
        dG0_std = []
        for pH in pH_range:
            dG0_prime, dG0_std, sqrt_Sigma = p['prediction_model'].get_transformed_dG0(ph=pH, ionic_strength=I, temperature=T)
            dG0_primes.append(dG0_prime)
            
        dG0_std = np.matrix(dG0_std)
        cc_dict['%s.dG0_primes' % entry] = dG0_primes
        cc_dict['%s.dG0_std' % entry] = dG0_std.T
        
    savemat('res/%s.mat' % file_name, cc_dict)
