import logging
from collections import Counter

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

logger = logging.getLogger(__name__)


class AbstractModel(dict):

    @classmethod
    def from_atoms(cls, atoms: Atoms, extra_info=None, **kwargs):
        """ASE's original implementation"""
        arrays = atoms.arrays.copy()
        n_atoms = len(atoms)
        dct = {
            'arrays': {
                'numbers': arrays.pop('numbers').tolist(),
                'positions': arrays.pop('positions').tolist(),

            },
            'info': {
                'cell': atoms.cell.tolist(),
                'pbc': atoms.pbc.tolist(),
                'constraints': [],
            },
        }

        for key, value in arrays.items():

            if isinstance(value, np.ndarray):
                dct['arrays'][key] = value.tolist()
                continue

            dct['arrays'][key] = value

        for key, value in atoms.info.items():

            if isinstance(value, np.ndarray):
                dct['info'][key] = value.tolist()
                continue

            dct['info'][key] = value

        if atoms.calc is not None:
            dct['info']['calculator_name'] = atoms.calc.__class__.__name__
            dct['info']['calculator_parameters'] = atoms.calc.todict()

            for key, value in atoms.calc.results.items():

                if isinstance(value, np.ndarray):
                    if value.shape[0] == n_atoms:
                        dct['arrays'][key] = value.tolist()
                    else:
                        dct['info'][key] = value.tolist()

                    continue

                dct['info'][key] = value

        # if atoms.constraints:
        #     dct['constraints'] = [c.todict() for c in atoms.constraints]

        if extra_info is not None:
            dct['info'].update(extra_info)

        dct['derived'] = {
            'elements': Counter(atoms.get_chemical_symbols()),
            'arrays_keys': list(dct['arrays'].keys()),
            'info_keys': list(dct['info'].keys())
        }

        return cls(**dct)

    def to_atoms(self):

        cell = self['info'].pop('cell', None)
        pbc = self['info'].pop('pbc', None)

        numbers = self['arrays'].pop('numbers', None)
        positions = self['arrays'].pop('positions', None)

        atoms = Atoms(numbers=numbers,
                      cell=cell,
                      pbc=pbc,
                      positions=positions)

        if 'calculator_name' in self['info']:
            calculator_name = self['info'].pop('calculator_name')
            params = self['info'].pop('calculator_parameters', {})
            results = self.pop('results', {})
            # TODO: Proper initialisation fo Calculators
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            atoms.calc = SinglePointCalculator(atoms, **params, **results)

        atoms.arrays.update(self['arrays'])
        atoms.info.update(self['info'])

        return atoms
