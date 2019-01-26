"""
Extended XYZ reader
"""

import re
from pathlib import Path
from typing import Union

re_key = r'(?P<key>[A-Za-z_]+[A-Za-z0-9_-]*)'
re_values = r'(?:' \
            r'(?P<single_value>[^\s"]+)' \
            r'|' \
            r'["\{\}](?P<quoted_value>[^"\{\}]+)["\{\}]' \
            r')\s*'

dtype_map = {
    'R': float,
    'I': int,
    'S': str,
    'L': lambda x: True if x in ['T', 'True'] else False
}


class XYZReader(object):
    REGEXP = re.compile(re_key + r's*=\s*' + re_values)

    def __init__(self, file: Union[str, Path]):
        self.file = file
        self._fp = None

    def __iter__(self):
        while self._fp is not None:
            yield self.read_frame()

    def __enter__(self):
        self._fp = open(self.file, 'r').__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self._fp.__exit__(exc_type, exc_val, exc_tb)

    def read_frame(self):
        n_atoms = self._read_number_of_atoms()
        frame_properties, property_types = self._read_frame_properties()
        atoms_properties = self._read_atoms_properties(n_atoms, property_types)

        return n_atoms, frame_properties, atoms_properties

    def _read_number_of_atoms(self):
        n_atoms = int(next(self._fp))
        return n_atoms

    def _read_frame_properties(self):
        line = next(self._fp)

        data = {key: single_value or quoted_value for key, single_value, quoted_value in self.REGEXP.findall(line)}
        property_types = data.pop('Properties')
        frame_properties = {key: convert(value) for key, value in data.items()}

        return frame_properties, property_types

    def _read_atoms_properties(self, n_atoms: int, property_types: str):
        props = property_types.split(':')  # Format is "[NAME:TYPE:NCOLS]...]", e.g. "species:S:1:pos:R:3".
        property_types = tuple((name, dtype_map[dtype], int(ncols))
                               for name, dtype, ncols
                               in zip(props[0::3], props[1::3], props[2::3]))

        results = {key: [] for key, _, _ in property_types}

        for _ in range(n_atoms):
            data = next(self._fp).split()
            ind = 0

            for (name, dtype, ncols) in property_types:
                if ncols == 1:
                    results[name].append(data[ind])
                else:
                    results[name].append([dtype(x) for x in data[ind:ind + ncols]])

                ind += ncols

        return results


def convert(text: str):
    """Convert string to python object by guessing its type"""

    # array-like object
    elements = text.split()
    if len(elements) > 1:
        return [convert(el) for el in elements]

    try:
        return int(text)
    except ValueError:
        pass

    try:
        return float(text)
    except ValueError:
        pass

    if text == 'T' or text == 'True':
        return True
    elif text == 'F' or text == 'False':
        return False
    else:
        # return as it is
        return text


if __name__ == '__main__':
    from pathlib import Path

    directory = Path('data/')
    file = directory / 'bcc_bulk_54_expanded_2_high.xyz'
    # file = directory / 'GAP_6.xyz'

    with XYZReader(file) as reader:
        for atoms in reader:
            n_atoms, frame_properties, atoms_properties = atoms
            # print('==========================')
            print(f'natoms: {n_atoms}')
            print(frame_properties)
            print(atoms_properties)
