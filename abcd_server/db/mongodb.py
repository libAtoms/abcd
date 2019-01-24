import json
import base64
import datetime
import numpy as np
from pprint import pprint

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from bson.objectid import ObjectId
from pymongo import MongoClient

from abcd_server.encoders.json import JSONEncoderOld, JSONDecoderOld

from abcd_server.encoders import JSONEncoder, DictEncoder
from abcd_server.db.base import Database

import ase
from os import linesep


class MongoDatabase(Database):
    """Wrapper to make database operations easy"""

    def __init__(self, url='mongodb://localhost:27017/', db_name='abcd', collection_name='atoms'):
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

    def push(self, atoms: ase.Atoms):
        pass

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
        return f'ABCD(type={self.__class__.__name__}, url={self.client.HOST}:{self.client.PORT}, ...)'

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


if __name__ == '__main__':
    from ase.io import iread

    database_name = 'fekad_test'
    collection_name = 'atoms'

    client = MongoClient('mongodb://fekad:qwe123@ds211613.mlab.com:11613/fekad_test')
    db = client[database_name]
    collection = db[collection_name]

    d = Database('mongodb://fekad:qwe123@ds211613.mlab.com:11613/fekad_test', db_name=database_name,
                 collection_name=collection_name)
    d.info()
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

    collection.insert_one(data)
