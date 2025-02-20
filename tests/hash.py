import io

import ase
from ase.calculators.lj import LennardJones
import ase.io
import numpy as np
import pytest

from abcd.model import AbstractModel


@pytest.fixture
def rng():
    return np.random.default_rng(seed=42)


def test_hash_structure(rng):
    # create atoms & add a calculator
    atoms = ase.Atoms(
        "H3",
        positions=rng.random(size=(3, 3)),
        pbc=True,
        cell=[2, 2, 2],
    )
    atoms.calc = LennardJones()
    atoms.calc.calculate(atoms)

    # dump to XYZ
    buffer = io.StringIO()
    ase.io.write(buffer, atoms, format="extxyz")

    # read back
    buffer.seek(0)
    atoms_read = ase.io.read(buffer, format="extxyz")
    assert atoms.positions[0, 0] == pytest.approx(atoms_read.positions[0, 0], abs=1e-10)

    # read in both of them
    abcd_data = AbstractModel.from_atoms(atoms)
    abcd_data_after_read = AbstractModel.from_atoms(atoms_read)

    assert atoms_read.positions[0, 0:5] == pytest.approx(
        atoms.positions[0, 0:5], abs=1e-10
    )
    assert abcd_data["hash_structure"] == abcd_data_after_read["hash_structure"], (
        atoms.positions,
        atoms_read.positions,
    )
