import pytest
from rdkit import Chem

from polygnn_kit import __version__
from polygnn_kit import polygnn_kit as pk


def test_version():
    assert __version__ == "0.1.0"


@pytest.fixture
def example_polymer_data():
    return pk.LinearPol("[*]=CC=[*]")


def test_periodic_bond_type(example_polymer_data):
    assert example_polymer_data.periodic_bond_type == Chem.BondType.DOUBLE


def test_two_star_bonds():
    """
    This test checks that creating a LinearPol from the following SMILES
    string will fail. This should happen because one star
    contains two bonds.
    """
    sm = "[*]CC/C=C\[*](Cl)"
    with pytest.raises(
        ValueError,
        match="Invalid repeat unit. It is likely that at least one star contains more than one bond.",
    ):
        pk.LinearPol(sm)


def test_mismatching_bonds():
    """
    This test checks that creating a LinearPol from the following SMILES
    string will fail. This should happen because one periodic bond is a
    single bond while the other periodic bond is a double bond.
    """
    sm = "[*]C1CCC(C=[*])C1"
    with pytest.raises(
        ValueError, match="Invalid repeat unit. Periodic bond types are mismatching."
    ):
        pk.LinearPol(sm)


@pytest.fixture
def example_ladder_data():
    sm = "C1C([g])C([e])OC([t])C1[d]"
    pol = pk.LadderPolymer(sm)
    return {
        "pol": pol,
    }


def test_LadderPolymer(example_ladder_data):
    pol = example_ladder_data["pol"]
    # indices below are assigned as "A" since
    # 2 < 4
    assert pol.starA1_ind == 2
    assert pol.connectorA1_ind == 1
    assert pol.starA2_ind == 9
    assert pol.connectorA2_ind == 8
    # indices below are assigned as "B" since
    # 4 < 2
    assert pol.starB1_ind == 4
    assert pol.connectorB1_ind == 3
    assert pol.starB2_ind == 7
    assert pol.connectorB2_ind == 6