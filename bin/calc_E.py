# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25th 2015

@author: flamholz

Modified on Thursday Jan 4th 2018

@author: joaocardoso
"""

import argparse
import csv

import numpy as np

from component_contribution.core.reaction import Reaction
from component_contribution.prediction_model.model import ComponentContributionModel
from component_contribution.thermodynamics.constants import STANDARD_T, F

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate reduction potentials for a number of reactions.')
    parser.add_argument('infile', type=argparse.FileType(), help='path to input file containing a list of reactions')
    parser.add_argument('outfile', type=argparse.FileType('w'), help='path to output file')
    parser.add_argument('--ionic_strength', default=0.2, type=float, help='ionic strength in molar units.')
    parser.add_argument('--pH_min', default=5, type=float, help='lowest pH to produce E0 for.')
    parser.add_argument('--pH_max', default=9, type=float, help='highest pH to produce E0 for.')
    parser.add_argument('--pH_step', default=0.05, type=float, help='pH increment.')

    args = parser.parse_args()

    ionic_strength = args.ionic_strength
    temperature = STANDARD_T

    cc = ComponentContributionModel.init()
    pHs = np.arange(args.pH_min, args.pH_max + args.pH_step,
                    args.pH_step)

    reactions_and_energies = []
    reader = csv.reader(args.infile)
    for row in reader:
        formula = row[0].strip()
        reaction = Reaction.from_equation(formula)
        reaction_atom_bag = reaction.atom_bag()
        n_e = reaction_atom_bag.pop('e-', 0)
        if len(reaction_atom_bag) != 0:
            raise ValueError('This is not a half-reaction (i.e. cannot be balanced by adding e-)')

        dG0_r, u_r = cc.get_reaction_dg(reaction)
        E0s = []
        for pH in pHs:
            ddG0_r = reaction.transform_ddG0(ph=pH, ionic_strength=ionic_strength, temperature=temperature)
            dG0_r_prime = dG0_r + ddG0_r
            E0_prime = 1000 * -dG0_r_prime / (n_e * F)  # mV
            E0s.append(E0_prime)

        reactions_and_energies.append((row, E0s))

    header = ['reaction', 'reaction_description', 'E\'m', 'Type', 'Source']
    pH_header = ['pH %.1f (mV)' % pH for pH in pHs]
    header += pH_header
    writer = csv.writer(args.outfile)
    writer.writerow(header)
    for rxn_data, pH_E0 in reactions_and_energies:
        energies_fmted = ['%.2f' % e for e in pH_E0]
        writer.writerow(rxn_data + energies_fmted)
