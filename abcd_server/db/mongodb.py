import json
import base64
import datetime
import numpy as np
from os import linesep

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from pymongo import MongoClient
from abcd_server.db.base import Database

from abcd_server.encoders import JSONEncoder, DictEncoder
from abcd_server.encoders.json import JSONEncoderOld, JSONDecoderOld
from abcd_server.encoders.dictionary import to_dict
from bson.objectid import ObjectId


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

    def push(self, atoms: Atoms):
        # with DictEncoder() as encoder:
        #     data = encoder.encode(atoms)
        data = atoms2dict(atoms)
        self.collection.insert_one(data)

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


all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy']


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


def atoms2dict(atoms):
    """ASE's original implementation"""
    dct = {
        'numbers': atoms.numbers.tolist(),
        'positions': atoms.positions.tolist(),
        # 'unique_id': '{}'.format(randint(16 ** 31, 16 ** 32 - 1))
    }
    if atoms.cell.any():
        dct['pbc'] = atoms.pbc.tolist()
        dct['cell'] = atoms.cell.tolist()
    if atoms.has('initial_magmoms'):
        dct['initial_magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.has('initial_charges'):
        dct['initial_charges'] = atoms.get_initial_charges()
    if atoms.has('masses'):
        dct['masses'] = atoms.get_masses()
    if atoms.has('tags'):
        dct['tags'] = atoms.get_tags()
    if atoms.has('momenta'):
        dct['momenta'] = atoms.get_momenta()
    if atoms.constraints:
        dct['constraints'] = [c.todict() for c in atoms.constraints]
    if atoms.calc is not None:
        dct['calculator'] = atoms.calc.name.lower()
        dct['calculator_parameters'] = atoms.calc.todict()
        if len(atoms.calc.check_state(atoms)) == 0:
            for prop in all_properties:
                try:
                    x = atoms.calc.get_property(prop, atoms, False)
                except PropertyNotImplementedError:
                    pass
                else:
                    if x is not None:
                        dct[prop] = x.tolist()
    return dct


if __name__ == '__main__':
    from ase.io import iread
    from pprint import pprint

    db = MongoDatabase('mongodb://localhost:27017/')
    db.info()

    for at in iread('../../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
        # print(at)
        at.calc.results['forces'] = at.arrays['force']
        # at.arrays['force'] = None

        json_data = json.dumps(at, cls=JSONEncoderOld)
        print(json_data)

        atom_dict = json.loads(json_data, cls=JSONDecoderOld)
        print(atom_dict)

        print(at == atom_dict)

    with JSONEncoder() as encoder:
        data = encoder.encode(at)

    print(data)

    with DictEncoder() as encoder:
        data = encoder.encode(at)

    pprint(data)

    from bson.json_util import dumps, _json_convert

    dumps(data)

    db.collection.insert_one(atoms2dict(at))
