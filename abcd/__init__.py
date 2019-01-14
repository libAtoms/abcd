import hashlib
import logging
import json
import base64
import numpy as np

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

# import ase

logger = logging.getLogger(__name__)


class ABCD(object):
    def __init__(self, url='http://localhost'):
        #  todo: authentication(user, token), https
        self.url = url

    def push(self, atoms):
        # todo: list of Atoms, metadata(user, project, tags)
        message = json.dumps(atoms)
        message_hash = hashlib.md5(message.encode('utf-8')).hexdigest()
        logger.info(message)
        return message_hash

    def pull(self, query=None, properties=None):
        atoms = json.loads(message)
        return atoms

    def query(self, query_string):
        pass


class Message(object):
    pass


def dumps(atoms):
    """Wrapper for generating json from Atoms object"""
    return json.dumps(atoms, cls=AtomsEncoder)


def loads(string):
    """Wrapper for creating Atoms object from json"""
    obj = json.loads(string, cls=AtomsDecoder)
    atoms = Atoms(
        numbers=obj['numbers'],
        pbc=obj['pbc'],
        cell=obj['cell'],
        positions=obj['positions'],
    )
    atoms.set_calculator(
        SinglePointCalculator(
            atoms,
            forces=obj['forces'],
            energy=obj['energy']
        )
    )

    return atoms


class Massage(dict):
    @staticmethod
    def b64_to_numpy(s):
        """b64_to_numpy: Helper function for converting base64 string into numpy array using buffers"""
        if s is None:
            return None
        return np.frombuffer(base64.b64decode(s))

    @classmethod
    def from_string(cls, string):
        obj = json.loads(string)

        obj['numbers'] = obj.get('numbers', None),
        obj['pbc'] = obj.get('pbc', None),
        obj['cell'] = cls.b64_to_numpy(obj['cell']).reshape(-1, 3),
        obj['positions'] = cls.b64_to_numpy(obj['positions']).reshape(-1, 3),
        obj['forces'] = np.frombuffer(base64.b64decode(obj['forces'])).reshape(-1, 3),
        obj['energy'] = np.frombuffer(base64.b64decode(obj['energy']))[0]
        return obj

    @classmethod
    def from_atoms(cls, atoms):
        return cls(
            numbers=atoms.numbers.tolist(),
            pbc=atoms.pbc.tolist(),
            cell=cls.numpy_to_b64(atoms.cell),
            positions=cls.numpy_to_b64(atoms.positions),
            forces=cls.numpy_to_b64(atoms.calc.results['forces']),
            energy=cls.numpy_to_b64(atoms.calc.results['energy']),
        )

    @staticmethod
    def numpy_to_b64(array):
        """numpy_to_b64: Helper function for converting numpy array into base64 string without loosing precision"""
        if array is None:
            return None
        return base64.b64encode(array).decode('UTF-8')


def numpy_to_b64(array):
    """numpy_to_b64: Helper function for converting numpy array into base64 string without loosing precision"""
    if array is None:
        return None
    return base64.b64encode(array).decode('UTF-8')


def atoms_2_dict(atoms):
    return {
        'numbers': atoms.numbers.tolist(),
        'pbc': atoms.pbc.tolist(),
        'cell': numpy_to_b64(atoms.cell),
        'positions': numpy_to_b64(atoms.positions),
        'forces': numpy_to_b64(atoms.calc.results['forces']),
        'energy': numpy_to_b64(atoms.calc.results['energy']),
    }


class AtomsEncoder(json.JSONEncoder):

    @staticmethod
    def numpy_to_b64(array):
        """numpy_to_b64: Helper function for converting numpy array into base64 string without loosing precision"""
        if array is None:
            return None
        return base64.b64encode(array).decode('UTF-8')

    def default(self, obj):
        if isinstance(obj, Atoms):
            # todo: masses, charges, formula
            return {
                # 'formula': obj.get_chemical_formula(),
                'numbers': obj.numbers.tolist(),
                'pbc': obj.pbc.tolist(),
                'cell': self.numpy_to_b64(obj.cell),
                'positions': self.numpy_to_b64(obj.positions),
                'forces': self.numpy_to_b64(obj.calc.results['forces']),
                'energy': self.numpy_to_b64(obj.calc.results['energy']),
            }
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


class AtomsDecoder(json.JSONDecoder):

    @staticmethod
    def b64_to_numpy(s):
        """b64_to_numpy: Helper function for converting base64 string into numpy array using buffers"""
        if s is None:
            return None
        return np.frombuffer(base64.b64decode(s))

    def decode(self, s, _w=None):
        obj = super().decode(s)
        obj['numbers'] = obj.get('numbers', None),
        obj['pbc'] = obj.get('pbc', None),
        obj['cell'] = self.b64_to_numpy(obj['cell']).reshape(-1, 3),
        obj['positions'] = self.b64_to_numpy(obj['positions']).reshape(-1, 3),
        obj['forces'] = np.frombuffer(base64.b64decode(obj['forces'])).reshape(-1, 3),
        obj['energy'] = np.frombuffer(base64.b64decode(obj['energy']))[0]
        return obj
        # atoms = Atoms(
        #     numbers=obj.get('numbers', None),
        #     pbc=obj.get('pbc', None),
        #     cell=self.b64_to_numpy(obj['cell']).reshape(-1, 3),
        #     positions=self.b64_to_numpy(obj['positions']).reshape(-1, 3),
        # )
        # atoms.set_calculator(
        #     SinglePointCalculator(
        #         atoms,
        #         forces=np.frombuffer(base64.b64decode(obj['forces'])).reshape(-1, 3),
        #         energy=np.frombuffer(base64.b64decode(obj['energy']))[0]
        #     )
        # )
        # return atoms


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    from ase.io import read

    # from ase.io import write
    # from ase.io.jsonio import read_json
    # write(atoms, '-', format='json')

    traj = read('../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None))

    for atoms in traj:
        # Hack to fix the representation of forces
        atoms.calc.results['forces'] = atoms.arrays['force']

    print(traj)

    atoms = traj[0]
    print(atoms)

    message = dumps(atoms)
    print(message)

    import requests
    payload = {'key1': 'value1', 'key2': 'value2'}
    r = requests.get('https://httpbin.org/get', params=payload)

    # new_atoms = loads(message)
    # print(new_atoms)
    #
    # assert atoms == new_atoms
    #
    # connection = ABCD()
    # for atoms in traj:
    #     hash_value = connection.push(atoms)
    #     print(hash_value)
