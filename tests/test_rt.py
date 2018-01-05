import logging
from component_contribution.core.compound import Compound

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)

for cid in ['C00002']:
    comp = Compound.from_kegg(cid)
    print(comp)
    print("%s: ddG0_prime = %.2f" % (cid, comp.transform(7.0, 0.2, 300)))
    print("%s: ddG0_prime(z=0) = %.2f" % (cid, comp.transform_neutral(7.0, 0.2, 300)))
