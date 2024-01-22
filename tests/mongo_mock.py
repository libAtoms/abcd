from io import StringIO
import logging
import unittest

from ase.io import read
from ase.atoms import Atoms
import mongomock

from abcd import ABCD


class MongoMock(unittest.TestCase):
    @classmethod
    @mongomock.patch(servers=(("localhost", 27017),))
    def setUpClass(cls):
        logging.basicConfig(level=logging.INFO)
        url = "mongodb://localhost"
        abcd = ABCD.from_url(url)
        abcd.print_info()

        cls.abcd = abcd

    @classmethod
    def tearDownClass(cls):
        cls.abcd.destroy()

    def test_thing(self):
        print(self.abcd.info())

    def test_push(self):
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

        self.abcd.destroy()
        self.abcd.push(atoms)
        new = list(self.abcd.get_atoms())[0]

        assert atoms == new
        self.abcd.destroy()


if __name__ == "__main__":
    unittest.main(verbosity=1, exit=False)
