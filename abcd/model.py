import datetime
import getpass
import logging
from collections import Counter, UserDict
from ase.calculators.singlepoint import SinglePointCalculator

import numpy as np
from ase import Atoms

logger = logging.getLogger(__name__)


class AbstractModel(UserDict):
    reserved_keys = {'n_atoms', 'cell', 'pbc', 'calculator_name', 'calculator_parameters'}

    def __init__(self, dict=None, **kwargs):
        self.arrays_keys = []
        self.info_keys = []
        self.results_keys = []
        self.derived_keys = []

        super().__init__(dict, **kwargs)

    @property
    def derived(self):
        return {
            'arrays_keys': self.arrays_keys,
            'info_keys': self.info_keys,
            'results_keys': self.results_keys,
            'derived_keys': self.derived_keys
        }

    def __getitem__(self, key):
        logger.info(f'__getitem__: {key}')

        if key == 'derived':
            return self.derived

        return super().__getitem__(key)

    def __setitem__(self, key, value):
        logger.info(f'__setitem__: {key}: {value}')

        if key == 'derived':
            raise KeyError('Please do not use "derived" as key because it is protected!')

        # if key not in self.derived:
        #     self.derived.append(key)

        super().__setitem__(key, value)

    def __delitem__(self, key):
        # self.derived.remove(key)
        super().__delitem__(key)

    def __iter__(self):
        for item in super().__iter__():
            yield item
        yield 'derived'

    @classmethod
    def from_atoms(cls, atoms: Atoms, extra_info=None, store_calc=True):
        """ASE's original implementation"""

        reserved_keys = {'n_atoms', 'cell', 'pbc', 'calculator_name', 'calculator_parameters', 'derived'}
        arrays_keys = set(atoms.arrays.keys())
        info_keys = set(atoms.info.keys())
        results_keys = set(atoms.calc.results.keys()) if store_calc and atoms.calc else {}

        all_keys = (reserved_keys, arrays_keys, info_keys, results_keys)
        if len(set.union(*all_keys)) != sum(map(len, all_keys)):
            print(all_keys)
            raise ValueError('All the keys must be unique!')

        item = cls()

        n_atoms = len(atoms)

        dct = {
            'n_atoms': n_atoms,
            'cell': atoms.cell.tolist(),
            'pbc': atoms.pbc.tolist(),
        }

        info_keys.update({'n_atoms', 'cell', 'pbc'})

        for key, value in atoms.arrays.items():
            if isinstance(value, np.ndarray):
                dct[key] = value.tolist()
            else:
                dct[key] = value

        for key, value in atoms.info.items():
            if isinstance(value, np.ndarray):
                dct[key] = value.tolist()
            else:
                dct[key] = value

        if store_calc and atoms.calc:
            dct['calculator_name'] = atoms.calc.__class__.__name__
            dct['calculator_parameters'] = atoms.calc.todict()
            info_keys.update({'calculator_name', 'calculator_parameters'})

            for key, value in atoms.calc.results.items():

                if isinstance(value, np.ndarray):
                    if value.shape[0] == n_atoms:
                        arrays_keys.update(key)
                    else:
                        info_keys.update(key)
                    dct[key] = value.tolist()

        item.arrays_keys = list(arrays_keys)
        item.info_keys = list(info_keys)
        item.results_keys = list(results_keys)

        item.update(dct)

        if extra_info:
            item.update(extra_info)
            item.info_keys.append(*extra_info.keys())

        item.pre_save()
        return item

    def to_atoms(self):
        arrays_keys = set(self.arrays_keys)
        info_keys = set(self.info_keys)

        cell = self.pop('cell', None)
        pbc = self.pop('pbc', None)
        numbers = self.pop('numbers', None)
        positions = self.pop('positions', None)
        results_keys = self.derived['results_keys']

        info_keys -= {'cell', 'pbc'}
        arrays_keys -= {'numbers', 'positions'}

        atoms = Atoms(
            cell=cell,
            pbc=pbc,
            numbers=numbers,
            positions=positions)

        if 'calculator_name' in self:
            # calculator_name = self['info'].pop('calculator_name')
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            params = self.pop('calculator_parameters', {})

            atoms.calc = SinglePointCalculator(atoms, **params)
            atoms.calc.results.update((key, self[key]) for key in results_keys)

        atoms.arrays.update((key, self[key]) for key in arrays_keys)
        atoms.arrays.update((key, self[key]) for key in info_keys)

        return atoms

    def pre_save(self):
        self.derived_keys = ['elements', 'username', 'uploaded', 'modified']

        cell = self['cell']

        if cell:
            volume = abs(np.linalg.det(cell))  # atoms.get_volume()
            self['volume'] = volume
            self.derived_keys.append('volume')

            virial = self.get('virial')
            if virial:
                # pressure P = -1/3 Tr(stress) = -1/3 Tr(virials/volume)
                self['pressure'] = -1 / 3 * np.trace(virial / volume)
                self.derived_keys.append('pressure')

        # 'elements': Counter(atoms.get_chemical_symbols()),
        self['elements'] = Counter(str(element) for element in self['numbers'])

        self['username'] = getpass.getuser()

        if not self.get('uploaded'):
            self['uploaded'] = datetime.datetime.utcnow()

        self['modified'] = datetime.datetime.utcnow()


if __name__ == '__main__':
    import io
    from pprint import pprint
    from ase.io import read

    logging.basicConfig(level=logging.INFO)
    # from ase.io import jsonio

    xyz = io.StringIO("""2
        Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
        Si       0.00000000       0.00000000       0.00000000 
        Si       0.00000000       0.00000000       0.00000000 
        """)

    atoms = read(xyz, format='xyz')
    atoms.set_cell([1, 1, 1])

    # print(atoms)
    # print(atoms.arrays)
    # print(atoms.info)

    # pprint(AbstractModel.from_atoms(atoms))

    # pprint(jsonio.encode(atoms.arrays))
    # pprint(jsonio.encode(atoms.info))
    # pprint(jsonio.encode(atoms.cell))
    #
    pprint(AbstractModel.from_atoms(atoms))

    model = AbstractModel.from_atoms(atoms)
    print(model.to_atoms())
