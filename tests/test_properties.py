import os

from pandas import DataFrame
import pytest

from abcd.backends.atoms_properties import Properties


class TestProperties:
    """Testing properties data reader"""

    @pytest.fixture(autouse=True)
    def property(self):
        """Load example data file."""
        class_path = os.path.normpath(os.path.abspath(__file__))
        data_file = os.path.dirname(class_path) + "/data/examples.csv"
        return Properties(data_file)

    def test_dataframe(self, property):
        """
        Test data correctly stored in pandas DataFrame.
        """
        assert isinstance(property.df, DataFrame)
        assert len(property.df) == 3

    def test_specify_units(self, property):
        """
        Test units can be specified manually, if they match existing fields.
        """
        input_units_1 = {"Integers": "items", "Floating": "seconds"}
        properties_1 = Properties(
            data_file=property.data_file,
            units=input_units_1,
        )
        assert properties_1.units == input_units_1

        input_units_2 = {"Fake": "m"}
        with pytest.raises(ValueError):
            properties_1 = Properties(
                data_file=property.data_file,
                units=input_units_2,
            )

    def test_infer_units(self, property):
        """
        Test units can be inferred from field names.
        """
        properties = Properties(
            data_file=property.data_file,
            infer_units=True,
        )
        expected_units = {"Comma units": "m", "Bracket units": "s"}
        expected_fields = [
            "Text",
            "Integers",
            "Floating",
            "Boolean",
            "Missing data",
            "Comma units",
            "Bracket units",
        ]
        assert properties.units == expected_units
        assert list(properties.df.columns.values) == expected_fields

    def test_struct_file(self, property):
        """
        Test structure file names can be inferred from a field.
        """
        struct_file_template = "test_{struct_name}_file.txt"
        struct_name_label = "Text"
        properties_1 = Properties(
            data_file=property.data_file,
            store_struct_file=True,
            struct_file_template=struct_file_template,
            struct_name_label=struct_name_label,
        )
        expected_struct_files = [
            "test_Some_file.txt",
            "test_test_file.txt",
            "test_data_file.txt",
        ]
        assert isinstance(properties_1.struct_files, list)
        for i, file in enumerate(expected_struct_files):
            assert properties_1.struct_files[i] == file

        invalid_template = "invalid_template"
        with pytest.raises(ValueError):
            Properties(
                data_file=property.data_file,
                store_struct_file=True,
                struct_file_template=invalid_template,
                struct_name_label=struct_name_label,
            )

        invalid_label = "label"
        with pytest.raises(ValueError):
            Properties(
                data_file=property.data_file,
                store_struct_file=True,
                struct_file_template=struct_file_template,
                struct_name_label=invalid_label,
            )

    def test_to_list(self, property):
        """
        Test dataframe can be converted into a list of properties.
        """
        assert len(property.to_list()) == 3
        assert isinstance(property.to_list(), list)
        assert isinstance(property.to_list()[0], dict)
        expected_property = {
            "Text": "Some",
            "Integers": 1,
            "Floating": 0.01,
            "Boolean": True,
            "Missing data": "Missing",
            "Comma units, m": 0,
            "Bracket units (s)": 0,
        }
        assert property.to_list()[0] == expected_property

    def test_missing_data(self, property):
        """
        Test missing data is not included in properties.
        """
        expected_property = {
            "Text": "test",
            "Integers": 2,
            "Floating": 0.1,
            "Boolean": False,
            "Comma units, m": 1,
            "Bracket units (s)": 1,
        }
        assert property.to_list()[1] == expected_property

    def test_to_list_units(self, property):
        """
        Test units are included in properties when converting to a list.
        """
        properties_1 = Properties(
            data_file=property.data_file,
            infer_units=True,
        )
        expected_units = {"Comma units": "m", "Bracket units": "s"}
        expected_property = {
            "Text": "Some",
            "Integers": 1,
            "Floating": 0.01,
            "Boolean": True,
            "Missing data": "Missing",
            "Comma units": 0,
            "Bracket units": 0,
            "units": expected_units,
        }
        assert properties_1.to_list()[0] == expected_property
