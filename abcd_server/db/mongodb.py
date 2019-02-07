import json
from os import linesep

import types
from typing import Union, Iterable

import ase
from pymongo import MongoClient

from abcd_server.db.base import Database
from abcd_server.formats import DictEncoder
from abcd_server.formats.myjson import JSONEncoderOld, JSONDecoderOld, JSONEncoder


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


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

    def push(self, atoms: Union[ase.Atoms, Iterable], extra_info=None):
        # with DictEncoder() as encoder:
        #     data = encoder.encode(atoms)
        if isinstance(atoms, ase.Atoms):
            data = DictEncoder().encode(atoms)
            if extra_info is not None:
                data['info'].update(extra_info)
            self.collection.insert_one(data)

        if isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):
            data = DictEncoder().encode_many(atoms, extra_info)
            self.collection.insert_many(data)

    def upload(self, file):
        data = ase.io.iread(file)
        self.push(data)

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        raise NotImplementedError

    def query(self, query_string):
        raise NotImplementedError

    def search(self, query_string: str):
        raise NotImplementedError

    def get_atoms(self, dbfilter):
        from abcd_server.formats.dictionary import DictDecoder
        with DictDecoder() as decoder:
            for atoms in self.db.atoms.find(dbfilter):
                yield decoder.decode(atoms)

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
        # {'_id': 'calculator_name', 'count': 56},
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
        return f'{self.__class__.__name__}(' \
            f'url={self.client.HOST}:{self.client.PORT}, ' \
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

    db.collection.insert_one(DictEncoder().encode(atoms))
