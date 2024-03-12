from io import StringIO
import logging
import os
from time import sleep
import unittest

from ase.atoms import Atoms
from ase.io import read
from opensearchpy.exceptions import ConnectionError

from abcd import ABCD
from abcd.backends.atoms_opensearch import AtomsModel, OpenSearchDatabase


class OpenSearch(unittest.TestCase):
    """
    Testing live OpenSearch database functions.
    """

    @classmethod
    def setUpClass(cls):
        """
        Set up OpenSearch database connection.
        """
        if os.getenv("GITHUB_ACTIONS") != "true":
            raise unittest.SkipTest("Only runs via GitHub Actions")
        cls.security_enabled = os.getenv("security_enabled") == "true"
        cls.port = int(os.environ["port"])
        cls.host = "localhost"
        if os.environ["opensearch-version"] == "latest":
            cls.credential = "admin:myStrongPassword123!"
        else:
            cls.credential = "admin:admin"

        logging.basicConfig(level=logging.INFO)

        url = f"opensearch://{cls.credential}@{cls.host}:{cls.port}"
        try:
            abcd = ABCD.from_url(
                url,
                index_name="test_index",
                use_ssl=cls.security_enabled,
            )
        except (ConnectionError, ConnectionResetError):
            sleep(10)
            abcd = ABCD.from_url(
                url,
                index_name="test_index",
                use_ssl=cls.security_enabled,
            )

        assert isinstance(abcd, OpenSearchDatabase)
        cls.abcd = abcd

    @classmethod
    def tearDownClass(cls):
        """
        Delete index from OpenSearch database.
        """
        cls.abcd.destroy()

    def push_data(self):
        """
        Helper function to upload an example xyz file to the database.
        """
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

    def test_info(self):
        """
        Test printing database info.
        """
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
        """
        Test destroying database index.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.abcd.refresh()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))

        self.abcd.destroy()
        self.assertFalse(self.abcd.client.indices.exists("test_index"))

    def test_create(self):
        """
        Test creating database index.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.abcd.refresh()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))
        self.assertFalse(self.abcd.client.indices.exists("fake_index"))

    def test_push(self):
        """
        Test pushing atoms objects to database individually.
        """
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

    def test_delete(self):
        """
        Test deleting all documents from database.
        """
        self.push_data()
        self.push_data()

        self.assertEqual(self.abcd.count(), 2)
        self.abcd.delete()
        self.assertTrue(self.abcd.client.indices.exists("test_index"))
        self.abcd.refresh()
        self.assertEqual(self.abcd.count(), 0)

    def test_bulk(self):
        """
        Test pushing atoms object to database together.
        """
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
        """
        Test counting the number of documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()
        self.push_data()
        self.assertEqual(self.abcd.count(), 2)

    def test_property(self):
        """
        Test getting values of a property from the database.
        """
        self.abcd.destroy()
        self.abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" prop_1="test_1" pbc="F F F"
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
            Properties=species:S:1:pos:R:3 s="sadf" prop_2="test_2" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        self.abcd.push(atoms_2)

        self.abcd.refresh()
        prop = self.abcd.property("prop_1")
        expected_prop = ["test_1"]
        self.assertEqual(prop, expected_prop)

    def test_properties(self):
        """
        Test getting all properties from the database.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()
        props = self.abcd.properties()
        expected_props = {
            "info": ["_vtk_test", "cell", "formula", "n_atoms", "pbc", "s", "volume"],
            "derived": [
                "elements",
                "hash",
                "hash_structure",
                "modified",
                "uploaded",
                "username",
                "volume",
            ],
            "arrays": ["numbers", "positions"],
        }
        self.assertEqual(props, expected_props)

    def test_count_property(self):
        """
        Test counting values of specified properties from the database.
        """
        self.abcd.destroy()
        self.abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" prop_1="1" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_1 = read(xyz_1, format="extxyz")
        assert isinstance(atoms_1, Atoms)
        atoms_1.set_cell([1, 1, 1])
        self.abcd.push(atoms_1)

        xyz_2 = StringIO(
            """1
            Properties=species:S:1:pos:R:3 s="sadf" prop_2="2" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        self.abcd.push(atoms_2)

        self.abcd.refresh()
        self.assertEqual(self.abcd.count_property("prop_1"), {1: 1})
        self.assertEqual(self.abcd.count_property("n_atoms"), {1: 1, 2: 1})
        self.assertEqual(self.abcd.count_property("volume"), {1.0: 2})

    def test_count_properties(self):
        """
        Test counting appearences of each property in documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" prop_1="test_1" pbc="F F F"
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
            Properties=species:S:1:pos:R:3 s="sadf" prop_2="test_2" pbc="F F F"
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
            "prop_1": {"count": 1, "category": "info", "dtype": "scalar(str)"},
            "prop_2": {"count": 1, "category": "info", "dtype": "scalar(str)"},
            "cell": {"count": 2, "category": "info", "dtype": "array(float)"},
            "elements": {"count": 2, "category": "derived", "dtype": "scalar(dict)"},
            "formula": {"count": 2, "category": "info", "dtype": "scalar(str)"},
            "hash": {"count": 2, "category": "derived", "dtype": "scalar(str)"},
            "hash_structure": {
                "count": 2,
                "category": "derived",
                "dtype": "scalar(str)",
            },
            "modified": {"count": 2, "category": "derived", "dtype": "scalar(str)"},
            "n_atoms": {"count": 2, "category": "info", "dtype": "scalar(int)"},
            "numbers": {"count": 2, "category": "arrays", "dtype": "vector(int, N)"},
            "pbc": {"count": 2, "category": "info", "dtype": "vector(bool)"},
            "positions": {
                "count": 2,
                "category": "arrays",
                "dtype": "array(float, N x 3)",
            },
            "s": {"count": 2, "category": "info", "dtype": "scalar(str)"},
            "uploaded": {"count": 2, "category": "derived", "dtype": "scalar(str)"},
            "username": {"count": 2, "category": "derived", "dtype": "scalar(str)"},
            "volume": {"count": 2, "category": "derived", "dtype": "scalar(float)"},
        }

        self.assertEqual(props, expected_counts)

    def test_add_property(self):
        """
        Test adding a property to documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()
        self.abcd.add_property({"TEST_PROPERTY": "TEST_VALUE"})

        self.abcd.refresh()
        data = self.abcd.client.search(index="test_index")
        self.assertEqual(
            data["hits"]["hits"][0]["_source"]["TEST_PROPERTY"], "TEST_VALUE"
        )
        self.assertIn(
            "TEST_PROPERTY", data["hits"]["hits"][0]["_source"]["derived"]["info_keys"]
        )

    def test_rename_property(self):
        """
        Test renaming a property for documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()
        self.abcd.add_property({"TEST_PROPERTY": "TEST_VALUE"})
        self.abcd.refresh()
        self.abcd.rename_property("TEST_PROPERTY", "NEW_PROPERTY")
        self.abcd.refresh()

        data = self.abcd.client.search(index="test_index")
        self.assertEqual(
            data["hits"]["hits"][0]["_source"]["NEW_PROPERTY"], "TEST_VALUE"
        )

    def test_delete_property(self):
        """
        Test deleting a property from documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()

        self.abcd.add_property({"TEST_PROPERTY": "TEST_VALUE"})
        self.abcd.refresh()
        data = self.abcd.client.search(index="test_index")
        self.assertEqual(
            data["hits"]["hits"][0]["_source"]["TEST_PROPERTY"], "TEST_VALUE"
        )

        self.abcd.delete_property("TEST_PROPERTY")
        self.abcd.refresh()
        data = self.abcd.client.search(index="test_index")
        with self.assertRaises(KeyError):
            data["hits"]["hits"][0]["_source"]["TEST_PROPERTY"]
        self.assertNotIn(
            "TEST_PROPERTY", data["hits"]["hits"][0]["_source"]["derived"]["info_keys"]
        )

    def test_get_items(self):
        """
        Test getting a dictionary of values from documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()

        expected_items = {
            "_id": None,
            "n_atoms": 2,
            "numbers": [14, 14],
            "_vtk_test": "t _ e s t",
            "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            "cell": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "pbc": [False, False, False],
            "volume": 1.0,
            "hash_structure": None,
            "s": "sadf",
            "elements": {"14": 2},
            "uploaded": None,
            "formula": "Si2",
            "modified": None,
            "derived": {
                "info_keys": [
                    "s",
                    "n_atoms",
                    "_vtk_test",
                    "cell",
                    "pbc",
                    "formula",
                    "volume",
                ],
                "derived_keys": [
                    "elements",
                    "username",
                    "uploaded",
                    "modified",
                    "volume",
                    "hash_structure",
                    "hash",
                ],
                "arrays_keys": ["numbers", "positions"],
                "results_keys": [],
            },
            "hash": None,
            "username": None,
        }

        self.abcd.refresh()
        items = list(self.abcd.get_items())[0]

        for key in expected_items:
            if key not in [
                "_id",
                "hash_structure",
                "uploaded",
                "modified",
                "hash",
                "username",
            ]:
                if isinstance(expected_items[key], dict):
                    for dict_key in expected_items[key]:
                        if isinstance(expected_items[key][dict_key], list):
                            self.assertEqual(
                                set(expected_items[key][dict_key]),
                                set(items[key][dict_key]),
                            )
                        else:
                            self.assertEqual(
                                expected_items[key][dict_key], items[key][dict_key]
                            )
                else:
                    self.assertEqual(expected_items[key], items[key])

    def test_get_atoms(self):
        """
        Test getting values from documents in the database as Atoms objects.
        """
        self.abcd.destroy()
        self.abcd.create()
        self.push_data()
        expected_atoms = Atoms(symbols="Si2", pbc=False, cell=[1.0, 1.0, 1.0])
        self.assertEqual(expected_atoms, list(self.abcd.get_atoms())[0])

    def test_query(self):
        """
        Test querying documents in the database.
        """
        self.abcd.destroy()
        self.abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 s="sadf" prop_1="test_1" pbc="F F F"
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
            Properties=species:S:1:pos:R:3 s="sadf" prop_2="test_2" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        self.abcd.push(atoms_2)
        self.abcd.refresh()

        query_dict = {"match": {"n_atoms": 2}}
        query_all = "volume: [0 TO 10]"
        query_1 = "prop_1: *"
        query_2 = "prop_2: *"
        self.assertEqual(self.abcd.count(query_dict), 2)
        self.assertEqual(self.abcd.count(query_all), 2)
        self.assertEqual(self.abcd.count(query_1), 1)
        self.assertEqual(self.abcd.count(query_2), 1)


if __name__ == "__main__":
    unittest.main(verbosity=1, exit=False)
