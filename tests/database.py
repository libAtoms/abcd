import unittest
import mongomock
from openmock import openmock

from abcd import ABCD
import logging


class Mongo(unittest.TestCase):
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
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

        xyz = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
<<<<<<< HEAD
            Si       0.00000000       0.00000000       0.00000000 
            Si       0.00000000       0.00000000       0.00000000 
=======
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
>>>>>>> c962bfe (Apply black formatting)
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


class OpenSearch(unittest.TestCase):
    @classmethod
    @openmock
    def setUpClass(cls):
        from abcd.backends.atoms_opensearch import OpenSearchDatabase

        logging.basicConfig(level=logging.INFO)
        url = "opensearch://admin:admin@localhost:9200"
        abcd = ABCD.from_url(url, index_name="test_index", analyse_schema=False)
        assert isinstance(abcd, OpenSearchDatabase)
        cls.abcd = abcd

    @classmethod
    def tearDownClass(cls):
        cls.abcd.destroy()

    def test_destroy(self):
        self.assertTrue(self.abcd.client.indices.exists("test_index"))
        self.abcd.destroy()
        self.assertFalse(self.abcd.client.indices.exists("test_index"))
        return

    def test_create(self):
        self.abcd.destroy()
        self.abcd.create()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))
        self.assertFalse(self.abcd.client.indices.exists("fake_index"))

    def test_push(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms
        from abcd.backends.atoms_opensearch import AtomsModel

        self.abcd.destroy()
        self.abcd.create()
        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )
        atoms_1 = read(xyz_1, format="extxyz")
        assert isinstance(atoms_1, Atoms)
        atoms_1.set_cell([1, 1, 1])
        self.abcd.push(atoms_1)

        xyz_2 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
            W       0.00000000       0.00000000       0.00000000
            W       0.00000000       0.00000000       0.00000000
            """
        )
        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])

        result = AtomsModel(
            None,
            None,
            self.abcd.client.search(index="test_index")["hits"]["hits"][0]["_source"],
        ).to_ase()
        assert atoms_1 == result
        assert atoms_2 != result

    def test_bulk(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms
        from abcd.backends.atoms_opensearch import AtomsModel

        self.abcd.destroy()
        self.abcd.create()
        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )
        atoms_1 = read(xyz_1, format="extxyz")
        assert isinstance(atoms_1, Atoms)
        atoms_1.set_cell([1, 1, 1])

        xyz_2 = StringIO(
            """1
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            """
        )
        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])

        atoms_list = []
        atoms_list.append(atoms_1)
        atoms_list.append(atoms_2)
        self.abcd.push(atoms_list)
        assert self.abcd.count() == 2

        result_1 = AtomsModel(
            None,
            None,
            self.abcd.client.search(index="test_index")["hits"]["hits"][0]["_source"],
        ).to_ase()
        result_2 = AtomsModel(
            None,
            None,
            self.abcd.client.search(index="test_index")["hits"]["hits"][1]["_source"],
        ).to_ase()
        assert atoms_1 == result_1
        assert atoms_2 == result_2

    def test_count(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

        self.abcd.destroy()
        self.abcd.create()
        xyz = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms = read(xyz, format="extxyz")
        assert isinstance(atoms, Atoms)
        atoms.set_cell([1, 1, 1])
        self.abcd.push(atoms)
        self.abcd.push(atoms)
        assert self.abcd.count() == 2


if __name__ == "__main__":
    unittest.main(verbosity=1, exit=False)
