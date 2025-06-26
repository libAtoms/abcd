# Copyright (c) 2025.
# Authors: Elliott Kasoar, Tamás K. Stenczel
# This program is distributed under the MIT License, see LICENSE.md.

from io import StringIO

from ase.io import read
import mongomock
import pytest

from abcd import ABCD


@pytest.fixture
@mongomock.patch(servers=(("localhost", 27017),))
def abcd_mongodb():
    url = "mongodb://localhost"
    abcd = ABCD.from_url(url)
    abcd.print_info()

    return abcd


def test_thing(abcd_mongodb):
    print(abcd_mongodb.info())


def test_push(abcd_mongodb):
    xyz = StringIO(
        """2
Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
Si       0.00000000       0.00000000       0.00000000
Si       0.00000000       0.00000000       0.00000000
"""
    )

    atoms = read(xyz, format="extxyz")
    atoms.set_cell([1, 1, 1])

    abcd_mongodb.destroy()
    abcd_mongodb.push(atoms)
    new = list(abcd_mongodb.get_atoms())[0]

    assert atoms == new
