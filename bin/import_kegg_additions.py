import os

from pandas import read_csv
from component_contribution.compound_cache import CompoundCache
from component_contribution.core.compound import Compound

current_directory = os.path.dirname(__file__)
package_directory = os.path.join(current_directory, "..")
data_directory = os.path.join(package_directory, "data")

additions = read_csv(os.path.join(data_directory, "kegg_additions.tsv"), sep="\t")

del current_directory
del package_directory
del data_directory

cache = CompoundCache.get_instance()
assert isinstance(cache, CompoundCache)

for _, row in additions.iterrows():
    compound = Compound.from_inchi(row.inchi, "KEGG", "C%05d" % row.cid)
    if compound.inchi_key in cache or compound.compound_id in cache:
        cache.replace(compound)
    else:
        cache.add(compound)

cache.dump()
