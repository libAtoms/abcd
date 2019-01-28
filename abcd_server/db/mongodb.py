import json
import base64
import datetime
import numpy as np
from os import linesep
import types
from typing import Union, Generator
from ase import Atoms
import ase
from ase.calculators.singlepoint import SinglePointCalculator

from pymongo import MongoClient
from abcd_server.db.base import Database

from abcd_server.encoders import JSONEncoder, DictEncoder
from abcd_server.encoders.json import JSONEncoderOld, JSONDecoderOld
from abcd_server.encoders.dictionary import to_dict
from bson.objectid import ObjectId


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


# all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
#                   'charges', 'magmom', 'magmoms', 'free_energy']

class MongoDatabase(Database):
    """Wrapper to make database operations easy"""

    def __init__(self, url='mongodb://localhost:27017/', db_name='abcd', collection_name='atoms'):
        super().__init__()

        self.client = MongoClient(url)
        self.db = self.client[db_name]
        self.collection = self.db[collection_name]

    def info(self):
        return {
            'host': self.client.HOST,
            'port': self.client.PORT,
            'db': self.db.name,
            'collection': self.collection.name,
            'number of confs': self.collection.count()
        }

    def destroy(self):
        self.collection.remove()

    def push(self, atoms: Union[ase.Atoms, Generator, list]):
        # with DictEncoder() as encoder:
        #     data = encoder.encode(atoms)
        if isinstance(atoms, ase.Atoms):
            data = atoms2dict(atoms)
            self.collection.insert_one(data)

        if isinstance(atoms, types.GeneratorType):
            raise NotImplementedError('Generators')

        if isinstance(atoms, list):
            raise NotImplementedError()
            # self.collection.insert_many()

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        raise NotImplementedError

    def query(self, query_string):
        raise NotImplementedError

    def search(self, query_string: str):
        raise NotImplementedError

    def get_atoms(self, id: str) -> Atoms:
        raise NotImplementedError

    def __repr__(self):
        return f'{self.__class__.__name__}(' \
            f'url={self.client.HOST}:{self.client.PORT}, ' \
            f'db={self.db.name}, ' \
            f'collection={self.collection.name})'

    def _repr_html_(self):
        """jupyter notebook representation"""
        return '<b>ABCD MongoDB database</b>'

    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join([
            '{:=^50}'.format(' ABCD MongoDB '),
            '{:>10}: {}'.format('type', 'mongodb'),
            linesep.join('{:>10}: {}'.format(k, v) for k, v in self.info().items())
        ])

        print(out)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


def atoms2dict(atoms: ase.Atoms) -> dict:
    """ASE's original implementation"""
    arrays = atoms.arrays.copy()

    dct = {
        'cell': atoms.cell.tolist(),
        'pbc': atoms.pbc.tolist(),
        'numbers': arrays.pop('numbers').tolist(),
        'positions': arrays.pop('positions').tolist(),
        'arrays': {},
        'info': {},
        'results': {},
        'constraints': [],
    }

    for key, value in arrays.items():

        if isinstance(value, np.ndarray):
            dct['arrays'][key] = value.tolist()
            continue

        dct[key] = value

    for key, value in atoms.info.items():

        if isinstance(value, np.ndarray):
            dct['info'][key] = value.tolist()
            continue

        dct['info'][key] = value

    if atoms.calc is not None:
        dct['results']['calculator_name'] = atoms.calc.name.lower(),
        dct['results']['calculator_parameters'] = atoms.calc.todict()

        for key, value in atoms.calc.results.items():

            if isinstance(value, np.ndarray):
                dct['results'][key] = value.tolist()
                continue

            dct['results'][key] = value

    # if atoms.constraints:
    #     dct['constraints'] = [c.todict() for c in atoms.constraints]

    return dct


if __name__ == '__main__':
    from ase.io import iread
    from pprint import pprint

    db = MongoDatabase('mongodb://localhost:27017/')
    db.info()

    for atoms in iread('../../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
        # print(at)
        atoms.calc.results['forces'] = atoms.arrays['force']
        # at.arrays['force'] = None

        json_data = json.dumps(atoms, cls=JSONEncoderOld)
        print(json_data)

        atom_dict = json.loads(json_data, cls=JSONDecoderOld)
        print(atom_dict)

        print(atoms == atom_dict)

    with JSONEncoder() as encoder:
        data = encoder.encode(atoms)

    print(data)

    with DictEncoder() as encoder:
        data = encoder.encode(atoms)

    pprint(data)

    from bson.json_util import dumps, _json_convert

    dumps(data)

    db.collection.insert_one(atoms2dict(atoms))
