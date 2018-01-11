import pytest
from component_contribution.compound_cache import CompoundCache
from component_contribution.core.compound import Compound

KEGG_COMPOUNDS = {"C00001", "C00002", "C00009"}
CHEBI_COMPOUNDS = {"CHEBI:15377", "CHEBI:15422", "CHEBI:18367"}


@pytest.fixture(scope="function")
def compound_cache(tmpdir):
    cache = str(tmpdir.mkdir(".cache").join("compound.cache"))
    CompoundCache._forget_class_instance_reference()
    ccache = CompoundCache.get_instance(cache_file_name=cache, bootstrap=False)
    assert len(ccache) == 0
    return ccache


@pytest.fixture(params=KEGG_COMPOUNDS)
def kegg_compound_id(request):
    return request.param


@pytest.fixture(params=CHEBI_COMPOUNDS)
def chebi_compound_id(request):
    return request.param


def test_add_kegg_compound(compound_cache, kegg_compound_id):
    assert isinstance(compound_cache, CompoundCache)
    cpd = compound_cache.get_compound(kegg_compound_id)

    assert isinstance(cpd, Compound)
    assert cpd.compound_id == kegg_compound_id
    assert cpd.database == "KEGG"
    assert kegg_compound_id in compound_cache.compound_id2inchi_key
    assert cpd.inchi_key in compound_cache.inchi_key2compound_id
    assert kegg_compound_id in compound_cache
    assert cpd.inchi_key in compound_cache
    assert cpd in compound_cache


def test_add_chebi_compound(compound_cache, chebi_compound_id):
    assert isinstance(compound_cache, CompoundCache)
    cpd = compound_cache.get_compound(chebi_compound_id)

    assert isinstance(cpd, Compound)
    assert cpd.compound_id == chebi_compound_id
    assert cpd.database == "CHEBI"
    assert chebi_compound_id in compound_cache.compound_id2inchi_key
    assert cpd.inchi_key in compound_cache.inchi_key2compound_id
    assert chebi_compound_id in compound_cache
    assert cpd.inchi_key in compound_cache
    assert cpd in compound_cache


def add_repeated_compound(compound_cache):
    cpd = compound_cache.add("C00009")

    assert isinstance(cpd, Compound)
    assert cpd.compound_id == "C00009"
    assert cpd.database == "KEGG"
    assert "C00001" in compound_cache.compound_id2inchi_key
    assert cpd.inchi_key in compound_cache.inchi_key2compound_id
    assert "C00001" in compound_cache
    assert cpd.inchi_key in compound_cache
    assert cpd in compound_cache

    cpd2 = compound_cache.add("CHEBI:15422")

    assert cpd2 == cpd
