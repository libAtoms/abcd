from io import StringIO
import logging
import os

from ase.atoms import Atoms
from ase.io import read
from openmock import openmock
import pytest

from abcd import ABCD
from abcd.backends.atoms_opensearch import AtomsModel, OpenSearchDatabase


class TestOpenSearchMock:
    """Testing mock OpenSearch database functions."""

    @pytest.fixture(autouse=True)
    @openmock
    def abcd(self):
        """Set up database connection."""

        if "port" in os.environ:
            port = int(os.environ["port"])
        else:
            port = 9200
        host = "localhost"

        logging.basicConfig(level=logging.INFO)

        url = f"opensearch://admin:admin@{host}:{port}"
        opensearch_abcd = ABCD.from_url(url, index_name="test_index", use_ssl=False)
        assert isinstance(opensearch_abcd, OpenSearchDatabase)
        return opensearch_abcd

    def test_destroy(self, abcd):
        """
        Test destroying database index.
        """
        assert abcd.client.indices.exists("test_index") is True
        abcd.destroy()
        assert abcd.client.indices.exists("test_index") is False

    def test_create(self, abcd):
        """
        Test creating database index.
        """
        abcd.destroy()
        abcd.create()
        assert abcd.client.indices.exists("test_index") is True
        abcd.client.indices.exists("fake_index") is False

    def test_push(self, abcd):
        """
        Test pushing atoms objects to database individually.
        """
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

        result = AtomsModel(
            None,
            None,
            abcd.client.search(index="test_index")["hits"]["hits"][0]["_source"],
        ).to_ase()
        assert atoms_1 == result
        assert atoms_2 != result

    def test_bulk(self, abcd):
        """
        Test pushing atoms object to database together.
        """
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
        """
        Test counting the number of documents in the database.
        """
        abcd.destroy()
        abcd.create()
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
        abcd.push(atoms)
        assert abcd.count() == 2
