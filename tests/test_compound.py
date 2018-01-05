import pytest
from component_contribution.core.compound import Compound

kegg_ids = {"C00001", "C00002", "C00009"}

smiles = {
    "C00001": '[OH2]',
    "C00002": 'C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc23)O1)O)O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O',
    "C00009": 'OP(O)(O)=O'
}

inchis = {
    "C00001": 'InChI=1S/H2O/h1H2',
    "C00002": 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)'
              '27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)',
    "C00009": 'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)'
}


formulas = {
    "C00001": 'H2O',
    "C00002": 'C10H16N5O13P3',
    "C00009": 'H3O4P'
}

pkas = {
    "C00001": [],
    "C00002": [12.6, 7.42, 5.0, 3.29, 2.53, 0.9],
    "C00009": [12.9, 6.95, 1.8]
}

#major_ms_at_ph7 = {1, [0; 0; 1; 0; 0; 0; 0], [0; 1; 0; 0]};
charges = {
    "C00001": [0],
    "C00002": [-5, -4, -3, -2, -1, 0, 1],
    "C00009": [-3, -2, -1, 0]
}

number_of_protons = {
    "C00001": [2],
    "C00002": [11, 12, 13, 14, 15, 16, 17],
    "C00009": [0, 1, 2, 3]
}


@pytest.fixture(params=kegg_ids)
def compound(request):
    kegg_id = request.param
    return Compound.from_kegg(kegg_id)


def test_compound(compound):
    assert compound.smiles_ph_7 == Compound.smiles2smiles(smiles[compound.compound_id])
    assert compound.p_kas == pkas[compound.compound_id]
    assert compound.inchi == inchis[compound.compound_id]
    assert compound.charges == charges[compound.compound_id]
    assert compound.number_of_protons == number_of_protons[compound.compound_id]