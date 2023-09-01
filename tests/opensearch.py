import unittest
from abcd import ABCD
import logging

class OpenSearch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        import os
        if os.getenv("GITHUB_ACTIONS") != "true":
            raise unittest.SkipTest("Only runs via GitHub Actions")
        
        cls.security_enabled = os.getenv("security_enabled") == "true"
        cls.port = int(os.environ["port"])
        cls.host = "localhost"

        from abcd.backends.atoms_opensearch import OpenSearchDatabase

        logging.basicConfig(level=logging.INFO)
        url = f"opensearch://admin:admin@{cls.host}:{cls.port}"
        abcd = ABCD.from_url(
            url,
            index_name="test_index",
            analyse_schema=False,
            use_ssl=cls.security_enabled
        )
        assert isinstance(abcd, OpenSearchDatabase)
        cls.abcd = abcd

    @classmethod
    def tearDownClass(cls):
        cls.abcd.destroy()

    def test_info(self):
        self.abcd.destroy()
        self.abcd.create()
        self.abcd.refresh()
        self.abcd.print_info()

        info = {
            "host": self.host,
            "port": self.port,
            "db": "abcd",
            "index": "test_index",
            "number of confs": 0,
            "type": "opensearch",
        }
        self.assertEqual(self.abcd.info(), info)

    def test_destroy(self):
        self.abcd.create()
        self.abcd.refresh()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))

        self.abcd.destroy()
        self.assertFalse(self.abcd.client.indices.exists("test_index"))
        return
    
    def test_delete(self):
        self.abcd.create()
        self.abcd.refresh()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))

        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

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
        self.abcd.push(atoms)
        self.abcd.push(atoms)

        self.abcd.refresh()
        self.assertEqual(self.abcd.count(), 2)

        self.abcd.delete()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))
        self.abcd.refresh()
        self.assertEqual(self.abcd.count(), 0)
        return

    def test_create(self):
        self.abcd.destroy()
        self.abcd.create()

        self.abcd.refresh()
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
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t _ e s t" pbc="F F F"
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
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t _ e s t" pbc="F F F"
            W       0.00000000       0.00000000       0.00000000
            W       0.00000000       0.00000000       0.00000000
            """
        )
        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])

        self.abcd.refresh()
        result = AtomsModel(
            None,
            None,
            self.abcd.client.search(index="test_index")["hits"]["hits"][0]["_source"],
        ).to_ase()
        self.assertEqual(atoms_1, result)
        self.assertNotEqual(atoms_2, result)

    def test_bulk(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms
        from abcd.backends.atoms_opensearch import AtomsModel

        self.abcd.destroy()
        self.abcd.create()
        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t _ e s t" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )
        atoms_1 = read(xyz_1, format="extxyz")
        assert isinstance(atoms_1, Atoms)
        atoms_1.set_cell([1, 1, 1])

        xyz_2 = StringIO(
            """1
            Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t _ e s t" pbc="F F F"
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

        self.abcd.refresh()
        self.assertEqual(self.abcd.count(), 2)
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
        self.assertEqual(atoms_1, result_1)
        self.assertEqual(atoms_2, result_2)

    def test_count(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

        self.abcd.destroy()
        self.abcd.create()
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
        self.abcd.push(atoms)
        self.abcd.push(atoms)

        self.abcd.refresh()
        self.assertEqual(self.abcd.count(), 2)

    def test_property(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

        self.abcd.destroy()
        self.abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" test_prop_1="test_prop_1" pbc="F F F"
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
            Properties=species:S:1:pos:R:3 s="sadf" test_prop_2="test_2" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        self.abcd.push(atoms_2)

        self.abcd.refresh()
        prop = self.abcd.property('test_prop_1')
        expected_prop = ['test_prop_1']
        self.assertEqual(prop, expected_prop)

    def test_properties(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

        self.abcd.destroy()
        self.abcd.create()

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
        self.abcd.push(atoms)

        self.abcd.refresh()
        props = self.abcd.properties()
        expected_props = {
            'info': ['_vtk_test', 'cell', 'formula', 'n_atoms', 'pbc', 's', 'volume'],
            'derived': [
                'elements',
                'hash',
                'hash_structure',
                'modified',
                'uploaded',
                'username',
                'volume'
            ],
            'arrays': ['numbers', 'positions']
        }
        self.assertEqual(props, expected_props)

    def test_count_properties(self):
        from io import StringIO
        from ase.io import read
        from ase.atoms import Atoms

        self.abcd.destroy()
        self.abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" test_prop_1="test_1" pbc="F F F"
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
            Properties=species:S:1:pos:R:3 s="sadf" test_prop_2="test_2" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        self.abcd.push(atoms_2)

        self.abcd.refresh()
        props = self.abcd.count_properties()
        expected_counts = {
            'test_prop_1': {'count': 1, 'category': 'info', 'dtype': 'scalar(str)'},
            'test_prop_2': {'count': 1, 'category': 'info', 'dtype': 'scalar(str)'},
            'cell': {'count': 2, 'category': 'info', 'dtype': 'array(float)'},
            'elements': {'count': 2, 'category': 'derived', 'dtype': 'scalar(dict)'},
            'formula': {'count': 2, 'category': 'info', 'dtype': 'scalar(str)'},
            'hash': {'count': 2, 'category': 'derived', 'dtype': 'scalar(str)'},
            'hash_structure': {'count': 2, 'category': 'derived', 'dtype': 'scalar(str)'},
            'modified': {'count': 2, 'category': 'derived', 'dtype': 'scalar(str)'},
            'n_atoms': {'count': 2, 'category': 'info', 'dtype': 'scalar(int)'},
            'numbers': {'count': 2, 'category': 'arrays', 'dtype': 'vector(int, N)'},
            'pbc': {'count': 2, 'category': 'info', 'dtype': 'vector(bool)'},
            'positions': {'count': 2, 'category': 'arrays', 'dtype': 'array(float, N x 3)'},
            's': {'count': 2, 'category': 'info', 'dtype': 'scalar(str)'},
            'uploaded': {'count': 2, 'category': 'derived', 'dtype': 'scalar(str)'},
            'username': {'count': 2, 'category': 'derived', 'dtype': 'scalar(str)'},
            'volume': {'count': 2, 'category': 'derived', 'dtype': 'scalar(float)'}
        }

        self.assertEqual(props, expected_counts)

if __name__ == "__main__":
    unittest.main(verbosity=1, exit=False)
