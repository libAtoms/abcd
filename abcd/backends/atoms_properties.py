from __future__ import annotations
import pandas as pd
import numpy as np
from typing import Union
from pathlib import Path

class Properties():
    """
    Wrapper to identify and manipulate properties to be passed
    as extra_info to the database.

    Attributes
    ----------
    csv_file: Union[str, Path]
        Name or path to csv file containing properties.
    store_struct_file: bool
        Whether to construct a filename for each structure.
    struct_file_template: str
        Template string for path to files containing structure.
    struct_name_label: str
        Field name in csv file containing values for `struct_name`.
    df: pd.Dataframe
        Dataframe containing loaded property data from csv file.
    units: Union[dict, None], optional
        Units.
    """
    def __init__(
        self,
        csv_file: Union[str, Path],
        store_struct_file: bool = False,
        struct_file_template: Union[str, None] = None,
        struct_name_label: Union[str, None] = None,
        units: Union[dict, None] = None
    ):
        """
        Initialises class.

        Parameters
        ----------
        csv_file: Union[str, Path]
            Path or filename of csv file containing properties to be loaded.
        store_struct_file: bool, optional
            If true, use struct_file_template and struct_name_label to
            construct filename for each structure. Default is `False`.
        struct_file_template: Union[str, None], optiona
            Template string for path to files containing structure.
            Required only if store_struct_file is True.
            Template must contain `{struct_name}`, to ensure a unique file
            for each structure. Default is `None`.
        struct_name_label: Union[str, None], optional
            Field name in csv file containing values for `struct_name`.
            Required only if store_struct_file is True. Default is `None`.
        units: Union[dict, None], optional
            Units for fields in csv file. If unspecified, _separate_units()
            is used to identify units in field names. Default is `None`.
        """
        self.csv_file = csv_file
        self.store_struct_file = store_struct_file
        if self.store_struct_file:
            if struct_file_template is not None:
                self.struct_file_template = struct_file_template
            else:
                raise ValueError((
                    "`struct_file_template` must be specified if "
                    "store_struct_file is True."
                ))
            if struct_name_label is not None:
                self.struct_name_label = struct_name_label
            else:
                raise ValueError(
                    "`struct_name_label` must be specified if store_struct_file is True."
                )

        self.df = pd.read_csv(self.csv_file)
        self.df.replace({np.nan: None}, inplace=True)

        if units is not None:
            for key in units.keys():
                if key not in self.df.columns.values:
                    raise ValueError((
                        f"Invalid field name: {key}. Keys in `units` must "
                        f"correspond to field names in the loaded data."
                    ))
            self.units = units
        else:
            self._separate_units()

    def _separate_units(self):
        """
        Parse field names to determine units.
        """
        columns = []
        self.units = {}
        for column in list(self.df.columns.values):
            if "," in column:
                column_name = column.split(",")[0].strip()
                self.units[column_name] = column.split(",")[1].strip()
            elif "(" in column:
                column_name = column.split("(")[0].strip()
                self.units[column_name] = column.split("(")[1].strip()[:-1]
            else:
                column_name = column

            columns.append(column)

        self.df.columns = columns

    def get_struct_file(self, struct_name: str) -> str:
        """
            Evaluate struct_file_template to determine structure filename
            for current structure.

            Parameters
            ----------
            struct_name: str
                Name of current structure.

            Returns
            -------
            Filename for the current structure.
        """
        if struct_name is None:
            raise ValueError("`struct_name` must be specified")
        if "{struct_name}" not in self.struct_file_template:
            raise ValueError((
                f"'struct_name' must be a variable in the template file: "
                f"{self.struct_file_template}"
            ))
        else:
            return eval(f"f'{self.struct_file_template}'")

    def to_list(self) -> list[dict]:
        """
        Convert dataframe into list of properties for each structure.

        Returns
        -------
        List of property dictionaries for each structure in the dataframe.
        """
        properties_list = []
        self.struct_files = []
        for i in range(len(self.df)):
            properties = self.df.iloc[i].to_dict()
            properties["units"] = self.units

            if self.store_struct_file:
                try:
                    struct_name = self.df.iloc[i][self.struct_name_label]
                except KeyError as e:
                    raise ValueError((
                        f"{self.struct_name_label} is not a valid column in "
                        f"the data loaded."
                    ))
                struct_file = self.get_struct_file(struct_name)
                properties["struct_file"] = struct_file

            properties_list.append(
                {key: value for key, value in properties.items() if value is not None}
            )
        return properties_list