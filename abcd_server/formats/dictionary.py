import numpy as np

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from abcd_server.formats.base import BaseEncoder, BaseDecoder


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


from collections import Counter


class DictEncoder(BaseEncoder):
    # default_properties = ['numbers', 'positions']

    def __init__(self):
        super().__init__()

    def encode(self, atoms: Atoms, extra_info=None) -> dict:
        """ASE's original implementation"""
        arrays = atoms.arrays.copy()
        natoms = len(atoms)
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
                    if value.shape[0] == natoms:
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
        return dct

    # def encode(self, atoms: Atoms) -> dict:
    #     """ASE's original implementation"""
    #     arrays = atoms.arrays.copy()
    #     natoms = len(atoms)
    #     dct = {
    #         'arrays': [
    #             {'name': 'numbers', 'value': arrays.pop('numbers').tolist()},
    #             {'name': 'positions', 'value': arrays.pop('positions').tolist()},
    #         ],
    #         'info': {
    #             'cell': atoms.cell.tolist(),
    #             'pbc': atoms.pbc.tolist(),
    #             'constraints': [],
    #         },
    #         'pbc': atoms.pbc.tolist(),
    #     }
    #
    #     for key, value in arrays.items():
    #
    #         if isinstance(value, np.ndarray):
    #             dct['arrays'].append({'name': key, 'value': value.tolist()})
    #             continue
    #
    #         dct['arrays'].append({'name': key, 'value': value})
    #
    #     for key, value in atoms.info.items():
    #
    #         if isinstance(value, np.ndarray):
    #             dct['info'][key] = value.tolist()
    #             continue
    #
    #         dct['info'][key] = value
    #
    #     if atoms.calc is not None:
    #         dct['info']['calculator_name'] = atoms.calc.__class__.__name__
    #         dct['info']['calculator_parameters'] = atoms.calc.todict()
    #
    #         for key, value in atoms.calc.results.items():
    #
    #             if isinstance(value, np.ndarray):
    #                 if value.shape[0] == natoms:
    #                     dct['arrays'].append({'name': key, 'value': value.tolist()})
    #                 else:
    #                     dct['info'][key] = value.tolist()
    #
    #                 continue
    #
    #             dct['info'][key] = value
    #
    #     return dct

    def encode_many(self, traj, extra_info=None):
        for atoms in traj:
            yield self.encode(atoms, extra_info)


class DictDecoder(BaseDecoder):

    def __init__(self):
        super().__init__()

    def decode(self, data: dict):
        cell = data['info'].pop('cell', None)
        pbc = data['info'].pop('pbc', None)

        numbers = data['arrays'].pop('numbers', None)
        positions = data['arrays'].pop('positions', None)

        atoms = Atoms(numbers=numbers,
                      cell=cell,
                      pbc=pbc,
                      positions=positions)

        if 'calculator_name' in data['info']:
            calculator_name = data['info'].pop('calculator_name')
            params = data['info'].pop('calculator_parameters', {})
            results = data.pop('results', {})
            # TODO: Proper initialisation fo Calculators
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            atoms.calc = SinglePointCalculator(atoms, **params, **results)

        atoms.arrays.update(data['arrays'])
        atoms.info.update(data['info'])

        return atoms


#
# class DictEncoder(BaseEncoder):
#     # default_properties = ['numbers', 'positions']
#
#     def __init__(self):
#         super().__init__()
#
#     def visit_atoms(self, atoms):
#         """walk through on the structure of an Atoms object """
#         data = {
#             'numbers': self.visit_numbers(atoms),
#             'cell': self.visit_cell(atoms),
#             'pbc': self.visit_pbc(atoms),
#             'positions': self.visit_positions(atoms),
#             'forces': self.visit_forces(atoms),
#             'energy': self.visit_energy(atoms),
#             'metadata': self.visit_metadata(atoms)
#         }
#         return data
#
#     def visit_metadata(self, atoms):
#         # metadata = {}
#         # for key, val in atoms.arrays.items():
#         #     if key not in self.default_properties:
#         #         metadata[key] = self.convert(val)
#         return {key: self.convert(val) for key, val in atoms.arrays.items() if key not in self.default_properties}
#
#     def visit_numbers(self, atoms):
#         return super().visit_numbers(atoms).tolist()
#
#     def visit_cell(self, atoms):
#         return super().visit_cell(atoms).tolist()
#
#     def visit_pbc(self, atoms):
#         return super().visit_pbc(atoms).tolist()
#
#     def visit_positions(self, atoms):
#         return super().visit_positions(atoms).tolist()
#
#     def visit_forces(self, atoms):
#         return super().visit_forces(atoms).tolist()
#
#     def convert(self, data):
#         if isinstance(data, np.ndarray):
#             return data.tolist()
#         if isinstance(data, dict):
#             return {key: self.convert(val) for key, val in data.items()}
#         return data


# TODO: unique id should be a hash value


if __name__ == '__main__':
    from pathlib import Path
    from ase.io import iread

    direcotry = Path('../../utils/data/')

    file = direcotry / 'bcc_bulk_54_expanded_2_high.xyz'

    for atoms in iread(file.as_posix(), index=slice(1)):
        print(atoms)

        # Fixing force label
        atoms.calc.results['forces'] = atoms.arrays.pop('force')

        with DictEncoder() as encoder:
            d = encoder.encode(atoms)

        print(d)

        with DictDecoder() as decoder:
            new_atoms = decoder.decode(d)

        print(new_atoms)
        print(new_atoms == atoms)
        #
