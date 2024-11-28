import pytest

from io import StringIO
from ase.io import read
import numpy as np

from abcd.model import AbstractModel


@pytest.fixture
def extxyz_file():
    return StringIO(
        """2
        Properties=species:S:1:pos:R:3:forces:R:3 energy=-1 pbc="F T F" info="test"
        Si  0.0  0.0  0.0  0.4  0.6  -0.4
        Si  0.0  0.0  0.0  -0.1  -0.5  -0.6
        """
    )


def test_from_atoms(extxyz_file):
    """Test extracting data from ASE Atoms object."""
    expected_forces = np.array([[0.4, 0.6, -0.4], [-0.1, -0.5, -0.6]])
    expected_stress = np.array([-1.0, -1.0, -1.0, -2.1, 2.0, 1.8])

    atoms = read(extxyz_file, format="extxyz")
    atoms.calc.results["stress"] = expected_stress
    data = AbstractModel.from_atoms(atoms)

    # Test info
    info_keys = {
        "pbc",
        "n_atoms",
        "cell",
        "formula",
        "calculator_name",
        "calculator_parameters",
        "info",
    }
    assert info_keys == set(data.info_keys)
    assert data["pbc"] == [False, True, False]
    assert data["n_atoms"] == 2
    assert len(data["cell"]) == 3
    assert all(arr == [0.0, 0.0, 0.0] for arr in data["cell"])
    assert data["formula"] == "Si2"
    assert data["info"] == "test"

    # Test arrays
    assert {"numbers", "positions"} == set(data.arrays_keys)

    # Test results
    assert {"energy", "stress", "forces"} == set(data.results_keys)
    assert data["energy"] == -1
    assert data["forces"] == pytest.approx(expected_forces)
    assert data["stress"] == pytest.approx(expected_stress)

    # Test derived
    derived_keys = {
        "elements",
        "username",
        "uploaded",
        "modified",
        "volume",
        "hash",
        "hash_structure",
    }
    assert derived_keys == set(data.derived_keys)


def test_from_atoms_no_calc(extxyz_file):
    """Test extracting data from ASE Atoms object without results."""
    expected_stress = np.array([-1.0, -1.0, -1.0, -2.1, 2.0, 1.8])

    atoms = read(extxyz_file, format="extxyz")
    atoms.calc.results["stress"] = expected_stress
    data = AbstractModel.from_atoms(atoms, store_calc=False)

    # Test info
    assert {"pbc", "n_atoms", "cell", "formula", "info"} == set(data.info_keys)
    assert data["pbc"] == [False, True, False]
    assert data["n_atoms"] == 2
    assert len(data["cell"]) == 3
    assert all(arr == [0.0, 0.0, 0.0] for arr in data["cell"])
    assert data["formula"] == "Si2"
    assert data["info"] == "test"

    # Test arrays
    assert {"numbers", "positions"} == set(data.arrays_keys)

    # Test results
    results_keys = {
        "energy",
        "forces",
        "stress",
        "calculator_name",
        "calculator_parameters",
    }
    assert all(key not in data for key in results_keys)

    # Test derived
    derived_keys = {
        "elements",
        "username",
        "uploaded",
        "modified",
        "volume",
        "hash",
        "hash_structure",
    }
    assert derived_keys == set(data.derived_keys)


def test_to_ase(extxyz_file):
    """Test returning data to ASE Atoms object with results."""
    atoms = read(extxyz_file, format="extxyz")
    data = AbstractModel.from_atoms(atoms, store_calc=True)

    new_atoms = data.to_ase()

    # Test info set
    assert new_atoms.cell == pytest.approx(atoms.cell)
    assert new_atoms.pbc == pytest.approx(atoms.pbc)
    assert new_atoms.positions == pytest.approx(atoms.positions)
    assert new_atoms.numbers == pytest.approx(atoms.numbers)

    assert new_atoms.info["n_atoms"] == len(atoms)
    assert new_atoms.info["formula"] == atoms.get_chemical_formula()

    assert new_atoms.calc.results["energy"] == pytest.approx(
        atoms.calc.results["energy"]
    )
    assert new_atoms.calc.results["forces"] == pytest.approx(
        atoms.calc.results["forces"]
    )


def test_to_ase_no_results(extxyz_file):
    """Test returning data to ASE Atoms object without results."""
    atoms = read(extxyz_file, format="extxyz")
    data = AbstractModel.from_atoms(atoms, store_calc=False)

    new_atoms = data.to_ase()

    # Test info set
    assert new_atoms.cell == pytest.approx(atoms.cell)
    assert new_atoms.pbc == pytest.approx(atoms.pbc)
    assert new_atoms.positions == pytest.approx(atoms.positions)
    assert new_atoms.numbers == pytest.approx(atoms.numbers)

    assert new_atoms.info["n_atoms"] == len(atoms)
    assert new_atoms.info["formula"] == atoms.get_chemical_formula()

    assert new_atoms.calc is None
