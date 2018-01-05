# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 21:00:31 2014

@author: eladn

Modified on Thursday Jan 4th 2018

@author: joaocardoso
"""
import sys
from component_contribution.compound_cache import CompoundCache

compound_id = sys.argv[1]
ccache = CompoundCache.get_instance()
sys.stderr.write('removing %s from cache ...\n' % compound_id)
ccache.remove(compound_id)
sys.stderr.write('recalculating SMILES and pKa values ...\n')
comp = ccache.get_compound(compound_id)
sys.stderr.write('writing new data to cache ...\n')
ccache.dump()

d = comp.to_json_dict()
sys.stderr.write(''.join(['%20s : %s\n' % (k, v) for (k, v) in d.iteritems()]))
