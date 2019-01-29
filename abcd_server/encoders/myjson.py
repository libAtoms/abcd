import json
import base64
from abcd_server.encoders.base import BaseEncoderOld
import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator


def b64_to_numpy(s):
    """b64_to_numpy: Helper function for converting base64 string into numpy array using buffers"""
    return np.frombuffer(base64.b64decode(s))


def numpy_to_b64(array):
    """numpy_to_b64: Helper function for converting numpy array into base64 string without loosing precision"""
    return base64.b64encode(array).decode('UTF-8')


class JSONNumpyEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return base64.b64encode(obj).decode('UTF-8')

        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class JSONEncoderOld(BaseEncoderOld):

    def visit_atoms(self, atoms):
        """walk through on the structure of an Atoms object """
        data = {
            'numbers': self.visit_numbers(atoms),
            'cell': self.visit_cell(atoms),
            'pbc': self.visit_pbc(atoms),
            'positions': self.visit_positions(atoms),
            'forces': self.visit_forces(atoms),
            'energy': self.visit_energy(atoms),
            'metadata': {}
        }
        # return data
        return json.dumps(data, cls=JSONNumpyEncoder)

    def visit_numbers(self, atoms):
        return super().visit_numbers(atoms).tolist()

    def visit_cell(self, atoms):
        return numpy_to_b64(super().visit_cell(atoms))

    def visit_pbc(self, atoms):
        return super().visit_pbc(atoms).tolist()

    def visit_positions(self, atoms):
        return numpy_to_b64(super().visit_positions(atoms))

    def visit_forces(self, atoms):
        return numpy_to_b64(super().visit_forces(atoms))

    def visit_energy(self, atoms):
        return float(super().visit_energy(atoms))


class JSONEncoderOld(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, Atoms):
            # todo: masses, charges, formula
            default_properties = ['numbers', 'positions', 'force']
            json_data = {
                # 'formula': obj.get_chemical_formula(),
                'numbers': obj.numbers.tolist(),
                'pbc': obj.pbc.tolist(),
                'cell': obj.cell,
                'positions': obj.positions,
                'forces': obj.get_forces(),
                'energy': float(obj.get_potential_energy()),
                'metadata': {}
            }

            for key, val in obj.arrays.items():
                if key not in default_properties:
                    json_data['metadata'][key] = val

            return json_data

        if isinstance(obj, np.ndarray):
            return numpy_to_b64(obj)

        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class JSONDecoderOld(json.JSONDecoder):

    @staticmethod
    def b64_to_numpy(s):
        """b64_to_numpy: Helper function for converting base64 string into numpy array using buffers"""
        if s is None:
            return None
        return np.frombuffer(base64.b64decode(s))

    def decode(self, s, _w=None):
        obj = super().decode(s)

        atoms = Atoms(
            numbers=obj.get('numbers', None),
            pbc=obj.get('pbc', None),
            cell=self.b64_to_numpy(obj['cell']).reshape(-1, 3),
            positions=self.b64_to_numpy(obj['positions']).reshape(-1, 3),
        )
        atoms.set_calculator(
            SinglePointCalculator(
                atoms,
                forces=np.frombuffer(base64.b64decode(obj['forces'])).reshape(-1, 3),
                energy=obj['energy']
            )
        )
        return atoms


if __name__ == '__main__':

    from ase.io import iread
    from pprint import pprint

    for atoms in iread('../../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
        # print(at)
        atoms.calc.results['forces'] = atoms.arrays['force']
        print(atoms)

        json_data = json.dumps(atoms, cls=JSONEncoderOld)
        print(json_data)

        atom_dict = json.loads(json_data, cls=JSONDecoderOld)
        print(atom_dict)

        print(atoms == atom_dict)

    with JSONEncoderOld() as encoder:
        data = encoder.encode(atoms)

    print(data)

# class Visitor(object):
#     def __init__(self, atoms, visitor):
#         self.atoms = atoms
#         self.visitor = visitor
#
#     def visit(self):
#         return self.visitor.visit_atoms(self.atoms)
#
# visitor = Visitor(at, DictEncoder())
# visitor = Visitor(at, JSONEncoder())
# pprint(visitor.visit())


# class DictConverter(object):
#     default_properties = ['numbers', 'positions']
#
#     def encode(self, atoms):
#         data = {
#             # 'formula': obj.get_chemical_formula(),
#             'numbers': atoms.numbers.tolist(),
#             'pbc': atoms.pbc.tolist(),
#             'cell': atoms.cell,
#             'positions': atoms.positions,
#             'forces': atoms.get_forces(),
#             'energy': atoms.get_potential_energy(),
#             'metadata': {}
#         }
#         for key, val in atoms.arrays.items():
#             if key not in self.default_properties:
#                 data['metadata'][key] = val
#
#         return data
#
#     def decode(self, data):
#         atoms = Atoms(
#             numbers=data.get('numbers', None),
#             pbc=data.get('pbc', None),
#             cell=data.get('cell', None),
#             positions=data.get('positions', None),
#         )
#         atoms.set_calculator(
#             SinglePointCalculator(
#                 atoms,
#                 forces=data.get('forces', None),
#                 energy=data.get('energy', None)
#             )
#         )
