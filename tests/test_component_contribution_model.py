import pytest

from component_contribution.core.reaction import Reaction
from component_contribution.prediction_model.model import ComponentContributionModel

target_reactions = {
    "1": {"C00031": -1, "C00469": 2, "C00011": 2},
    "2": {"C00031": 1, "C00469": -2, "C00011": -2},
    "3": {"C00007": -1, "C00090": -2, "C02351": 2, "C00001": 1},
    "4": {"C00007": 1, "C00090": 2, "C02351": -2, "C00001": -1}
}

reactions_dg0 = {
    "1": -265.0,
    "2": 265.0,
    "3": -253.4,
    "4": 253.4
}


@pytest.fixture(scope='module')
def component_contribution_model():
    return ComponentContributionModel.init(matfile=None)


@pytest.fixture(scope='function', params=list(target_reactions))
def reaction(request):
    return Reaction(stoichiometry=target_reactions[request.param], database="eQuilibrator", reaction_id=request.param)


def test_component_contribution_model_prediction(component_contribution_model, reaction):
    dg0 = reactions_dg0[reaction.reaction_id]
    predicted, _ = component_contribution_model.get_reaction_dg(reaction)
    assert predicted == pytest.approx(dg0, 0.1)
