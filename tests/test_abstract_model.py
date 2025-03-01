# Copyright (c) 2025.
# Authors: Elliott Kasoar
# This program is distributed under the MIT License, see LICENSE.md.

import datetime
import io
from io import StringIO

import ase
from ase.calculators.lj import LennardJones
from ase.io import read, write
import numpy as np
import pytest
from pytest import approx

from abcd.model import AbstractModel, Hasher


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


def test_from_atoms_len_atoms_3():
    atoms = ase.Atoms(
        "H3",
        positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        pbc=True,
        cell=[2, 2, 2],
    )
    atoms.calc = LennardJones()
    atoms.calc.calculate(atoms)

    # convert
    abcd_data = AbstractModel.from_atoms(atoms, store_calc=True)

    assert set(abcd_data.info_keys) == {
        "pbc",
        "n_atoms",
        "cell",
        "formula",
        "calculator_name",
        "calculator_parameters",
    }
    assert set(abcd_data.arrays_keys) == {"numbers", "positions"}
    assert set(abcd_data.results_keys) == {
        "stress",
        "energy",
        "forces",
        "energies",
        "stresses",
        "free_energy",
    }

    # check some values as well
    assert abcd_data["energy"] == atoms.get_potential_energy()
    assert abcd_data["forces"] == approx(atoms.get_forces())


@pytest.mark.parametrize("store_calc", [True, False])
def test_write_and_read(store_calc):
    # create atoms & add a calculator
    atoms = ase.Atoms(
        "H3",
        positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        pbc=True,
        cell=[2, 2, 2],
    )
    atoms.calc = LennardJones()
    atoms.calc.calculate(atoms)

    # dump to XYZ
    buffer = io.StringIO()
    write(buffer, atoms, format="extxyz")

    # read back
    buffer.seek(0)
    atoms_read = read(buffer, format="extxyz")

    # read in both of them
    abcd_data = AbstractModel.from_atoms(atoms, store_calc=store_calc)
    abcd_data_after_read = AbstractModel.from_atoms(atoms_read, store_calc=store_calc)

    # check that all results are the same
    for key in ["info_keys", "arrays_keys", "derived_keys", "results_keys"]:
        assert set(getattr(abcd_data, key)) == set(
            getattr(abcd_data_after_read, key)
        ), f"{key} mismatched"

    # info & arrays same, except calc recognised as LJ when not from XYZ
    for key in set(abcd_data.info_keys + abcd_data.arrays_keys) - {
        "calculator_name",
        "calculator_parameters",
    }:
        assert abcd_data[key] == abcd_data_after_read[key], (
            f"{key}'s value does not match"
        )

    # date & hashed will differ
    for key in set(abcd_data.derived_keys) - {
        "hash",
        "modified",
        "uploaded",
    }:
        assert abcd_data[key] == abcd_data_after_read[key], (
            f"{key}'s value does not match"
        )

    # expected differences - n.b. order of calls above
    assert abcd_data_after_read["modified"] > abcd_data["modified"]
    assert abcd_data_after_read["uploaded"] > abcd_data["uploaded"]

    # expect results to match within fp precision
    for key in set(abcd_data.results_keys):
        assert abcd_data[key] == approx(np.array(abcd_data_after_read[key])), (
            f"{key}'s value does not match"
        )


def test_hash_update():
    """Test hash can be updated after initialisation."""
    hasher_1 = Hasher()

    init_hash = hasher_1()
    hasher_1.update("Test value")
    assert hasher_1() != init_hash


@pytest.mark.parametrize(
    "data",
    [
        1296,
        3.14,
        [1, 2, 3],
        (4, 5, 6),
        {"a": "value"},
        datetime.datetime.now(datetime.timezone.utc),
        b"test",
    ],
)
def test_hash_data_types(data):
    """Test updating hash for different data types."""
    hasher_1 = Hasher()
    hasher_1.update("Test value")
    updated_hash = hasher_1()

    hasher_1.update(data)
    assert updated_hash != hasher_1()


def test_second_hash_init():
    """Test second hash is initialised correctly."""
    hasher_1 = Hasher()

    init_hash = hasher_1()
    hasher_1.update("Test value")

    hasher_2 = Hasher()
    assert hasher_2() == init_hash


@pytest.mark.parametrize(
    "data",
    [
        1296,
        3.14,
        [1, 2, 3],
        (4, 5, 6),
        {"a": "value"},
        datetime.datetime.now(datetime.timezone.utc),
        b"test",
    ],
)
def test_consistent_hash(data):
    """Test two hashers agree with same data."""
    hasher_1 = Hasher()
    hasher_1.update("Test value")
    hasher_1.update(data)

    hasher_2 = Hasher()
    hasher_2.update("Test value")
    hasher_2.update(data)
    assert hasher_1() == hasher_2()
