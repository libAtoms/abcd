import types
from typing import Union, Iterable
from os import linesep

from pymongo import MongoClient

from ase import Atoms
from ase.io import iread

from abcd import ABCD
from abcd.model import AtomsModel
from abcd.parsers.arguments import key_val_str_to_dict


class MongoDatabase(ABCD):
    """Wrapper to make database operations easy"""

    def __init__(self, host='localhost', port=27017,
                 db_name='abcd', collection_name='atoms',
                 username=None, password=None, authSource='admin', **kwargs):
        super().__init__()

        print(host, port, db_name, collection_name, username, password, authSource, kwargs)

        self.client = MongoClient(
            host=host, port=port, username=username, password=password,
            authSource=authSource)

        self.db = self.client[db_name]
        self.collection = self.db[collection_name]

    def info(self):
        host, port = self.client.address

        return {
            'host': host,
            'port': port,
            'db': self.db.name,
            'collection': self.collection.name,
            'number of confs': self.collection.count(),
            'type': 'mongodb (pymongo)'
        }

    def destroy(self):
        self.collection.drop()

    def push(self, atoms: Union[Atoms, Iterable], extra_info=None):

        if extra_info and isinstance(extra_info, str):
            extra_info = key_val_str_to_dict(extra_info)

        if isinstance(atoms, Atoms):
            data = AtomsModel.from_atoms(atoms, extra_info)
            self.collection.insert_one(data)

        if isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):
            self.collection.insert_many(AtomsModel.from_atoms(at, extra_info) for at in atoms)

    def upload(self, file, extra_info=None):
        if extra_info:
            extra_info = key_val_str_to_dict(' '.join(extra_info))
        else:
            extra_info = {}

        extra_info['filename'] = str(file)
        data = iread(str(file))
        self.push(data, extra_info)

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

        return '{}('.format(self.__class__.__name__) + \
               'url={}:{}, '.format(host, port) + \
               'db={}, '.format(self.db.name) + \
               'collection={})'.format(self.collection.name)

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
    db = MongoDatabase(username='mongoadmin', password='secret')
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
