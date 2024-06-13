from io import StringIO
import logging
import os
from time import sleep

from ase.atoms import Atoms
from ase.io import read
from opensearchpy.exceptions import ConnectionError
import pytest

from abcd import ABCD
from abcd.backends.atoms_opensearch import AtomsModel, OpenSearchDatabase

NOT_GTHUB_ACTIONS = True
if os.getenv("GITHUB_ACTIONS") == "true":
    NOT_GTHUB_ACTIONS = False


@pytest.mark.skipif(NOT_GTHUB_ACTIONS, reason="Not running via GitHub Actions")
class TestOpenSearch:
    """Testing live OpenSearch database functions."""

    @pytest.fixture(autouse=True)
    def abcd(self):
        """Set up OpenSearch database connection."""
        security_enabled = os.getenv("security_enabled") == "true"
        self.port = int(os.environ["port"])
        self.host = "localhost"
        if os.environ["opensearch-version"] == "latest":
            credential = "admin:myStrongPassword123!"
        else:
            credential = "admin:admin"

        logging.basicConfig(level=logging.INFO)

        url = f"opensearch://{credential}@{self.host}:{self.port}"
        try:
            abcd_opensearch = ABCD.from_url(
                url,
                index_name="test_index",
                use_ssl=security_enabled,
            )
        except (ConnectionError, ConnectionResetError):
            sleep(10)
            abcd_opensearch = ABCD.from_url(
                url,
                index_name="test_index",
                use_ssl=security_enabled,
            )

        assert isinstance(abcd_opensearch, OpenSearchDatabase)
        return abcd_opensearch

    def push_data(self, abcd):
        """Helper function to upload an example xyz file to the database."""
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
        abcd.push(atoms)
        abcd.refresh()

    def test_info(self, abcd):
        """Test printing database info."""
        abcd.destroy()
        abcd.create()
        abcd.refresh()
        abcd.print_info()

        info = {
            "host": self.host,
            "port": self.port,
            "db": "abcd",
            "index": "test_index",
            "number of confs": 0,
            "type": "opensearch",
        }
        assert abcd.info() == info

    def test_destroy(self, abcd):
        """Test destroying database index."""
        abcd.destroy()
        abcd.create()
        abcd.refresh()
        assert abcd.client.indices.exists("test_index") is True

        abcd.destroy()
        assert abcd.client.indices.exists("test_index") is False

    def test_create(self, abcd):
        """Test creating database index."""
        abcd.destroy()
        abcd.create()
        abcd.refresh()
        assert abcd.client.indices.exists("test_index") is True
        assert abcd.client.indices.exists("fake_index") is False

    def test_push(self, abcd):
        """Test pushing atoms objects to database individually."""
        abcd.destroy()
        abcd.create()
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
        abcd.push(atoms_1)

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

        abcd.refresh()
        result = AtomsModel(
            None,
            None,
            abcd.client.search(index="test_index")["hits"]["hits"][0]["_source"],
        ).to_ase()
        assert atoms_1 == result
        assert atoms_2 != result

    def test_delete(self, abcd):
        """Test deleting all documents from database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)
        self.push_data(abcd)
        abcd.refresh()

        assert abcd.count() == 2
        abcd.delete()
        assert abcd.client.indices.exists("test_index") is True
        abcd.refresh()
        assert abcd.count() == 0

    def test_bulk(self, abcd):
        """Test pushing atoms object to database together."""
        abcd.destroy()
        abcd.create()
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
        abcd.push(atoms_list)

        abcd.refresh()
        assert abcd.count() == 2
        result_1 = AtomsModel(
            None,
            None,
            abcd.client.search(index="test_index")["hits"]["hits"][0]["_source"],
        ).to_ase()
        result_2 = AtomsModel(
            None,
            None,
            abcd.client.search(index="test_index")["hits"]["hits"][1]["_source"],
        ).to_ase()
        assert atoms_1 == result_1
        assert atoms_2 == result_2

    def test_count(self, abcd):
        """Test counting the number of documents in the database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)
        self.push_data(abcd)
        assert abcd.count() == 2

    def test_property(self, abcd):
        """Test getting values of a property from the database."""
        abcd.destroy()
        abcd.create()

        xyz_1 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 energy=-5.0 prop_1="test_1"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_1 = read(xyz_1, format="extxyz")
        assert isinstance(atoms_1, Atoms)
        atoms_1.set_cell([1, 1, 1])
        abcd.push(atoms_1, store_calc=False)

        xyz_2 = StringIO(
            """2
            Properties=species:S:1:pos:R:3 energy=-10.0 prop_2="test_2"
            Si       0.00000000       0.00000000       0.00000000
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        abcd.push(atoms_2, store_calc=False)

        abcd.refresh()
        prop = abcd.property("prop_1")
        expected_prop = ["test_1"]
        assert prop == expected_prop

        prop = abcd.property("energy")
        expected_prop = [-5.0, -10.0]
        assert prop[0] == expected_prop[0]
        assert prop[1] == expected_prop[1]

    def test_properties(self, abcd):
        """Test getting all properties from the database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)
        props = abcd.properties()
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
        assert props == expected_props

    def test_count_property(self, abcd):
        """Test counting values of specified properties from the database."""
        abcd.destroy()
        abcd.create()

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
        abcd.push(atoms_1)

        xyz_2 = StringIO(
            """1
            Properties=species:S:1:pos:R:3 s="sadf" prop_2="2" pbc="F F F"
            Si       0.00000000       0.00000000       0.00000000
            """
        )

        atoms_2 = read(xyz_2, format="extxyz")
        assert isinstance(atoms_2, Atoms)
        atoms_2.set_cell([1, 1, 1])
        abcd.push(atoms_2)

        abcd.refresh()
        assert abcd.count_property("prop_1") == {1: 1}
        assert abcd.count_property("n_atoms") == {1: 1, 2: 1}
        assert abcd.count_property("volume") == {1.0: 2}

    def test_count_properties(self, abcd):
        """Test counting appearences of each property in documents in the database."""
        abcd.destroy()
        abcd.create()

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
        abcd.push(atoms_1)

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
        abcd.push(atoms_2)

        abcd.refresh()
        props = abcd.count_properties()
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

        assert props == expected_counts

    def test_add_property(self, abcd):
        """Test adding a property to documents in the database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)
        abcd.add_property({"TEST_PROPERTY": "TEST_VALUE"})

        abcd.refresh()
        data = abcd.client.search(index="test_index")
        assert data["hits"]["hits"][0]["_source"]["TEST_PROPERTY"] == "TEST_VALUE"
        assert (
            "TEST_PROPERTY"
            in data["hits"]["hits"][0]["_source"]["derived"]["info_keys"]
        )

    def test_rename_property(self, abcd):
        """Test renaming a property for documents in the database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)
        abcd.add_property({"TEST_PROPERTY": "TEST_VALUE"})
        abcd.refresh()
        abcd.rename_property("TEST_PROPERTY", "NEW_PROPERTY")
        abcd.refresh()

        data = abcd.client.search(index="test_index")
        assert data["hits"]["hits"][0]["_source"]["NEW_PROPERTY"] == "TEST_VALUE"

    def test_delete_property(self, abcd):
        """Test deleting a property from documents in the database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)

        abcd.add_property({"TEST_PROPERTY": "TEST_VALUE"})
        abcd.refresh()
        data = abcd.client.search(index="test_index")
        assert data["hits"]["hits"][0]["_source"]["TEST_PROPERTY"] == "TEST_VALUE"

        abcd.delete_property("TEST_PROPERTY")
        abcd.refresh()
        data = abcd.client.search(index="test_index")
        with pytest.raises(KeyError):
            data["hits"]["hits"][0]["_source"]["TEST_PROPERTY"]
        assert (
            "TEST_PROPERTY"
            not in data["hits"]["hits"][0]["_source"]["derived"]["info_keys"]
        )

    def test_get_items(self, abcd):
        """Test getting a dictionary of values from documents in the database."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)

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

        abcd.refresh()
        items = list(abcd.get_items())[0]

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
                            assert set(expected_items[key][dict_key]) == set(
                                items[key][dict_key]
                            )
                        else:
                            assert expected_items[key][dict_key] == items[key][dict_key]
                else:
                    assert expected_items[key] == items[key]

    def test_get_atoms(self, abcd):
        """Test getting values from documents in the database as Atoms objects."""
        abcd.destroy()
        abcd.create()
        self.push_data(abcd)
        expected_atoms = Atoms(symbols="Si2", pbc=False, cell=[1.0, 1.0, 1.0])
        assert expected_atoms == list(abcd.get_atoms())[0]

    def test_query(self, abcd):
        """Test querying documents in the database."""
        abcd.destroy()
        abcd.create()

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
        abcd.push(atoms_1)

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
        abcd.push(atoms_2)
        abcd.refresh()

        query_dict = {"match": {"n_atoms": 2}}
        query_all = "volume: [0 TO 10]"
        query_1 = "prop_1: *"
        query_2 = "prop_2: *"
        assert abcd.count(query_dict) == 2
        assert abcd.count(query_all) == 2
        assert abcd.count(query_1) == 1
        assert abcd.count(query_2) == 1
