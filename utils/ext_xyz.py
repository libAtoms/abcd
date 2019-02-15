"""
Extended XYZ reader
"""

import re
from pathlib import Path
from typing import Union
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

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

        property_types = data.pop('Properties', None)
        # Default set of properties is atomic symbols and positions only
        if property_types is None:
            property_types = 'species:S:1:pos:R:3'

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
                    ind += ncols
                    continue

                results[name].append([dtype(x) for x in data[ind:ind + ncols]])
                ind += ncols

        return results

    def read_atoms(self, energy_label='energy', forces_label='forces'):
        for n_atoms, frame_properties, atoms_properties in self:

            cell = frame_properties.pop('Lattice', None)
            if cell is not None:
                cell = np.array(cell, order='F').reshape((3, 3))

            pbc = frame_properties.pop('pbc', None)
            if cell is not None and pbc is None:
                pbc = [True, True, True]

            symbols = atoms_properties.pop('species', None)
            numbers = atoms_properties.pop('Z', None)
            if symbols is not None:
                symbols = [s.capitalize() for s in symbols]
                numbers = None

            positions = atoms_properties.pop('pos', None)
            charges = atoms_properties.pop('charge', None)

            atoms = Atoms(symbols=symbols,
                          positions=positions,
                          numbers=numbers,
                          charges=charges,
                          cell=cell,
                          pbc=pbc)

            # Load results of previous calculations into SinglePointCalculator
            results = {}

            energy = frame_properties.pop(energy_label, None)
            if energy is not None:
                results['energy'] = energy

            magmom = frame_properties.pop('magmom', None)
            if magmom is not None:
                results['magmom'] = magmom

            free_energy = frame_properties.pop('free_energy', None)
            if free_energy is not None:
                results['free_energy'] = free_energy

            forces = atoms_properties.pop(forces_label, None)
            if forces is not None:
                results['forces'] = forces

            dipole = atoms_properties.pop('dipole', None)
            # TODO: Make sure that it has the proper representation
            if dipole is not None:
                results['dipole'] = dipole

            charges = atoms_properties.pop('charges', None)
            # TODO: Make sure that it has the proper representation
            if charges is not None:
                results['charges'] = charges

            magmoms = atoms_properties.pop('magmoms', None)
            # TODO: Make sure that it has the proper representation
            if magmoms is not None:
                results['magmoms'] = magmoms

            stress = atoms_properties.pop('stress', None)
            if stress is not None:
                stress = np.array(stress).reshape((3, 3), order='F')
                stress = np.array([stress[0, 0],
                                   stress[1, 1],
                                   stress[2, 2],
                                   stress[1, 2],
                                   stress[0, 2],
                                   stress[0, 1]])
                results['stress'] = stress

            if results:
                calculator = SinglePointCalculator(atoms, **results)
                atoms.set_calculator(calculator)

            # Storing all the remaining properties in the info
            atoms.info.update(frame_properties)
            atoms.info.update(atoms_properties)

            yield atoms


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

    directory = Path('../tutorials/data/')
    file = directory / 'bcc_bulk_54_expanded_2_high.xyz'
    # file = directory / 'GAP_6.xyz'

    with XYZReader(file) as reader:
        for atoms in reader.read_atoms(forces_label='force'):
            print('==========================')
            print(atoms)

    with XYZReader(file) as reader:
        for atoms in reader:
            n_atoms, frame_properties, atoms_properties = atoms
            # print('==========================')
            print(f'natoms: {n_atoms}')
            print(frame_properties)
            print(atoms_properties)
