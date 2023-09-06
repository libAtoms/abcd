import unittest
from abcd.backends.atoms_properties import Properties


class PropertiesTests(unittest.TestCase):
    """Testing properties data reader"""

    @classmethod
    def setUpClass(cls):
        """
        Load example data file.
        """
        import os

        class_path = os.path.normpath(os.path.abspath(__file__))
        data_file = os.path.dirname(class_path) + "/examples.csv"
        cls.property = Properties(data_file)

    def test_dataframe(self):
        """
        Test data correctly stored in pandas DataFrame.
        """
        from pandas import DataFrame

        assert isinstance(self.property.df, DataFrame)
        assert len(self.property.df) == 3

    def test_specify_units(self):
        """
        Test units can be specified manually, if they match existing fields.
        """
        input_units_1 = {"Integers": "items", "Floating": "seconds"}
        properties_1 = Properties(
            data_file=self.property.data_file,
            units=input_units_1,
        )
        self.assertEqual(properties_1.units, input_units_1)

        input_units_2 = {"Fake": "m"}
        with self.assertRaises(ValueError):
            properties_1 = Properties(
                data_file=self.property.data_file,
                units=input_units_2,
            )

    def test_infer_units(self):
        """
        Test units can be inferred from field names.
        """
        properties = Properties(
            data_file=self.property.data_file,
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
        self.assertEqual(properties.units, expected_units)
        self.assertEqual(list(properties.df.columns.values), expected_fields)

    def test_struct_file(self):
        """
        Test structure file names can be inferred from a field.
        """
        struct_file_template = "test_{struct_name}_file.txt"
        struct_name_label = "Text"
        properties_1 = Properties(
            data_file=self.property.data_file,
            store_struct_file=True,
            struct_file_template=struct_file_template,
            struct_name_label=struct_name_label,
        )
        expected_struct_files = [
            "test_Some_file.txt",
            "test_test_file.txt",
            "test_data_file.txt",
        ]
        self.assertIsInstance(properties_1.struct_files, list)
        for i, file in enumerate(expected_struct_files):
            self.assertEqual(properties_1.struct_files[i], file)

        invalid_template = "invalid_template"
        with self.assertRaises(ValueError):
            Properties(
                data_file=self.property.data_file,
                store_struct_file=True,
                struct_file_template=invalid_template,
                struct_name_label=struct_name_label,
            )

        invalid_label = "label"
        with self.assertRaises(ValueError):
            Properties(
                data_file=self.property.data_file,
                store_struct_file=True,
                struct_file_template=struct_file_template,
                struct_name_label=invalid_label,
            )

    def test_to_list(self):
        """
        Test dataframe can be converted into a list of properties.
        """
        self.assertEqual(len(self.property.to_list()), 3)
        self.assertIsInstance(self.property.to_list(), list)
        self.assertIsInstance(self.property.to_list()[0], dict)
        expected_property = {
            "Text": "Some",
            "Integers": 1,
            "Floating": 0.01,
            "Boolean": True,
            "Missing data": "Missing",
            "Comma units, m": 0,
            "Bracket units (s)": 0,
        }
        self.assertEqual(self.property.to_list()[0], expected_property)

    def test_missing_data(self):
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
        self.assertEqual(self.property.to_list()[1], expected_property)

    def test_to_list_units(self):
        """
        Test units are included in properties when converting to a list.
        """
        properties_1 = Properties(
            data_file=self.property.data_file,
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
        self.assertEqual(properties_1.to_list()[0], expected_property)


if __name__ == "__main__":
    unittest.main()
