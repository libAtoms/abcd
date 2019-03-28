from os import linesep

import types
from typing import Union, Iterable
from collections import Counter

from ase import Atoms
from ase.io import iread
from ase.calculators.singlepoint import SinglePointCalculator

import numpy as np

from pymongo import MongoClient

from abcd.backends.abstract import Database


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


class AtomsModel(dict):

    @classmethod
    def from_atoms(cls, atoms: Atoms, extra_info=None, **kwargs):
        dct = cls()

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

        return dct

    def to_atoms(self):
        data = self

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


class MongoDatabase(Database):
    """Wrapper to make database operations easy"""

    def __init__(self, url='mongodb://localhost:27017/', db='abcd', collection='atoms', **kwargs):
        super().__init__()

        self.client = MongoClient(url, **kwargs)
        self.db = self.client[db]
        self.collection = self.db[collection]

    def info(self):
        host, port = self.client.address

        return {
            'host': host,
            'port': port,
            'db': self.db.name,
            'collection': self.collection.name,
            'number of confs': self.collection.count()
        }

    def destroy(self):
        self.collection.remove()

    def push(self, atoms: Union[Atoms, Iterable], extra_info=None):

        if isinstance(atoms, Atoms):
            data = AtomsModel.from_atoms(atoms)
            if extra_info is not None:
                data['info'].update(extra_info)
            self.collection.insert_one(data)

        if isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):
            self.collection.insert_many(AtomsModel.from_atoms(at) for at in atoms)

    def upload(self, file):
        data = iread(file)
        self.push(data)

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        raise NotImplementedError

    def query(self, query_string):
        raise NotImplementedError

    def search(self, query_string: str):
        raise NotImplementedError

    def get_atoms(self, dbfilter):
        for dct in self.db.atoms.find(dbfilter):
            yield AtomsModel(dct).to_atoms()

    def count(self, dbfilter=None):
        if dbfilter is None:
            return self.collection.count()

        return self.db.atoms.count_documents(dbfilter)

    def get_property(self, name, dbfilter=None):

        if dbfilter is None:
            dbfilter = {}

        pipeline = [
            {'$match': dbfilter},
            {'$match': {'{}'.format(name): {"$exists": True}}},
            {'$project': {'_id': False, 'data': '${}'.format(name)}}
        ]

        return [val['data'] for val in self.collection.aggregate(pipeline)]

    def count_properties(self, dbfilter=None):

        if dbfilter is None:
            dbfilter = {}

        properties = {
            'info': {},
            'arrays': {}
        }

        pipeline = [
            {'$match': dbfilter},
            {'$unwind': '$derived.info_keys'},
            {'$group': {'_id': '$derived.info_keys', 'count': {'$sum': 1}}}
        ]

        info_keys = self.db.atoms.aggregate(pipeline)

        for val in info_keys:
            properties['info'][val['_id']] = {'count': val['count']}

        pipeline = [
            {'$match': dbfilter},
            {'$unwind': '$derived.arrays_keys'},
            {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}
        ]
        arrays_keys = list(self.db.atoms.aggregate(pipeline))
        for val in arrays_keys:
            properties['arrays'][val['_id']] = {'count': val['count']}

        return properties

    def __repr__(self):
        host, port = self.client.address

        return f'{self.__class__.__name__}(' \
            f'url={host}:{port}, ' \
            f'db={self.db.name}, ' \
            f'collection={self.collection.name})'

    def _repr_html_(self):
        """Jupyter notebook representation"""
        return '<b>ABCD MongoDB database</b>'

    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join(['{:=^50}'.format(' ABCD MongoDB '),
                            '{:>10}: {}'.format('type', 'mongodb'),
                            linesep.join('{:>10}: {}'.format(k, v) for k, v in self.info().items())])

        print(out)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


if __name__ == '__main__':
    # import json
    # from ase.io import iread
    # from pprint import pprint
    # from server.styles.myjson import JSONEncoderOld, JSONDecoderOld, JSONEncoder

    print('hello')
    db = MongoDatabase('mongodb://localhost:27017/')
    print(db.info())

    # for atoms in iread('../../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
    #     # print(at)
    #     atoms.calc.results['forces'] = atoms.arrays['force']
    #     # at.arrays['force'] = None
    #
    #     json_data = json.dumps(atoms, cls=JSONEncoderOld)
    #     print(json_data)
    #
    #     atom_dict = json.loads(json_data, cls=JSONDecoderOld)
    #     print(atom_dict)
    #
    #     print(atoms == atom_dict)
    #
    # with JSONEncoder() as encoder:
    #     data = encoder.encode(atoms)
    #
    # print(data)
    #
    # with DictEncoder() as encoder:
    #     data = encoder.encode(atoms)
    #
    # pprint(data)
    #
    # db.collection.insert_one(DictEncoder().encode(atoms))
