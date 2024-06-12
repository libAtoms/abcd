from io import StringIO
import logging
import unittest

from ase.io import read
from ase.atoms import Atoms
import mongomock
import pytest

from abcd import ABCD


class TestMongoMock:
    @pytest.fixture(autouse=True)
    @mongomock.patch(servers=(("localhost", 27017),))
    def abcd(self):
        logging.basicConfig(level=logging.INFO)
        url = "mongodb://localhost"
        mongo_abcd = ABCD.from_url(url)
        mongo_abcd.print_info()
        return mongo_abcd

    def test_info(self, abcd):
        print(abcd.info())

    def test_push(self, abcd):
        xyz = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t _ e s t" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms = read(xyz, format="extxyz")
        assert isinstance(atoms, Atoms)
        atoms.set_cell([1, 1, 1])

        abcd.destroy()
        abcd.push(atoms)
        new = list(abcd.get_atoms())[0]

        assert atoms == new
        abcd.destroy()
