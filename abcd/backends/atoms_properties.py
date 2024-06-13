from __future__ import annotations
import pandas as pd
import numpy as np
from typing import Union
from pathlib import Path
import chardet


class Properties:
    """
    Wrapper to identify and manipulate properties to be passed
    as extra_info to the database.

    Attributes
    ----------
    data_file: Union[str, Path]
        Name or path to data file containing properties. Treated as a csv file
        by default, but Excel spreadsheets may also be read.
    store_struct_file: bool
        Whether to construct a filename for each structure.
    struct_file_template: str
        Template string for path to files containing structure.
    struct_name_label: str
        Field name in data file containing values for `struct_name`.
    df: pd.Dataframe
        Dataframe containing loaded property data from data file.
    units: Union[dict, None], optional
        Units.
    struct_files: list[str]
        List containing a filename for each structure in the dataframe.
    encoding: str, optional
        Encoding of csv file to be read. Default is `utf-8`.
    """

    def __init__(
        self,
        data_file: Union[str, Path],
        store_struct_file: bool = False,
        struct_file_template: Union[str, None] = None,
        struct_name_label: Union[str, None] = None,
        units: Union[dict, None] = None,
        infer_units: bool = False,
        encoding: str = "utf-8",
    ):
        """
        Initialises class.

        Parameters
        ----------
        data_file: Union[str, Path]
            Path or filename of data file containing properties to be loaded.
            Assumed to be a csv file by default, but Excel spreadsheets may
            also be read.
        store_struct_file: bool, optional
            If true, use struct_file_template and struct_name_label to
            construct filename for each structure. Default is `False`.
        struct_file_template: Union[str, None], optiona
            Template string for path to files containing structure.
            Required only if store_struct_file is True.
            Template must contain `{struct_name}`, to ensure a unique file
            for each structure. Default is `None`.
        struct_name_label: Union[str, None], optional
            Field name in data file containing values for `struct_name`.
            Required only if store_struct_file is True. Default is `None`.
        units: Union[dict, None], optional
            Units for fields in data file. If unspecified, _separate_units()
            is used to identify units in field names. Default is `None`.
        infer_units: bool, optional
            Whether to attempt to infer units from field names in the
            dataframe. Unused if units is not `None`. Default is `False`.
        encoding: str, optional
            Encoding of file to be read. Default is `utf-8`.
            For pandas==1.2, setting this to `None` means `errors='replace'`
            is passed to `open()`, which replaces invalid characters with
            the replacement character. Otherwise, `errors='strict'` is passed
            to `open()`, which means UnicodeDecodeError are thrown if the
            encoding is wrong.
            For pandas==1.3, `encoding` no longer defines how errors are
            handled. `encoding_errors` instead defaults to `strict`, which has
            the same effect as non-None values of `encoding` for pandas==1.2.
        """
        self.data_file = data_file
        self.encoding = encoding
        try:
            self.df = pd.read_csv(self.data_file, encoding=self.encoding)
        except UnicodeDecodeError:
            detected = chardet.detect(Path(self.data_file).read_bytes())
            raise ValueError(
                f"File cannot be decoded using encoding: {self.encoding}."
                f" Detected encoding: {detected}."
            )
        except pd.errors.ParserError:
            self.df = pd.read_excel(self.data_file, header=0)

        self.df.replace({np.nan: None}, inplace=True)

        if units is not None:
            for key in units:
                if key not in self.df.columns.values:
                    raise ValueError(
                        f"Invalid field name: {key}. Keys in `units` must "
                        "correspond to field names in the loaded data."
                    )
            self.units = units
        elif infer_units:
            self._separate_units()
        else:
            self.units = None

        self.store_struct_file = store_struct_file
        if self.store_struct_file:
            if struct_file_template is None:
                raise ValueError(
                    "`struct_file_template` must be specified if "
                    "store_struct_file is True."
                )
            self.struct_file_template = struct_file_template

            if struct_name_label is None:
                raise ValueError(
                    "`struct_name_label` must be specified if store_struct_file is"
                    " True."
                )
            self.struct_name_label = struct_name_label
            self.set_struct_files()

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

            columns.append(column_name)

        self.df.columns = columns

    def set_struct_files(self):
        """
        Sets a list containing a filename for each structure in the dataframe.
        """
        self.struct_files = []

        for i in range(len(self.df)):
            try:
                struct_name = self.df.iloc[i][self.struct_name_label]
            except KeyError:
                raise ValueError(
                    f"{self.struct_name_label} is not a valid column in "
                    "the data loaded."
                )
            struct_file = self.get_struct_file(struct_name)
            self.struct_files.append(struct_file)

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
        if "{struct_name}" not in self.struct_file_template:
            raise ValueError(
                "'struct_name' must be a variable in the template file: "
                f"{self.struct_file_template}"
            )
        return eval(f"f'{self.struct_file_template}'")

    def to_list(self) -> list[dict]:
        """
        Convert dataframe into list of properties for each structure.

        Returns
        -------
        List of property dictionaries for each structure in the dataframe.
        """
        properties_list = []
        for i in range(len(self.df)):
            properties = self.df.iloc[i].to_dict()
            if self.units is not None:
                properties["units"] = self.units
            properties_list.append(
                {key: value for key, value in properties.items() if value is not None}
            )
        return properties_list
