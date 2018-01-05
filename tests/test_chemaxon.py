import pytest
from collections import namedtuple
from component_contribution.chemaxon import formula_and_charge, dissociation_constants


Chemical = namedtuple("Chemical", ['name', 'inchi'])

molecules = {
    "D-Erythrulose": 'InChI=1S/C4H8O4/c5-1-3(7)4(8)2-6/h3,5-7H,1-2H2/t3-/m1/s1'
}

npsas = {
    "D-Erythrulose": [-8.53, -3.86, -3.33, -3.01, 12.54, 13.9, 15.45]
}

charges = {
    "D-Erythrulose": 0
}

formulas = {
    "D-Erythrulose": "C4H8O4"
}


@pytest.fixture(scope='module', params=list(molecules.keys()))
def chemical(request):
    name = request.param
    inchi = molecules[name]
    return Chemical(name, inchi)


def test_formula_and_charge(chemical: Chemical):
    formula, charge = formula_and_charge(chemical.inchi)
    assert formula == formulas[chemical.name]
    assert charge == charges[chemical.name]


def test_dissociation_constants(chemical: Chemical):
    dissociation_table, major_ms = dissociation_constants(chemical.inchi)
    assert dissociation_table == npsas[chemical.name]
