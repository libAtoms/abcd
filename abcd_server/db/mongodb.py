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


class Database(object):
    def __init__(self, url='mongodb://localhost:27017/', db_name='abcd', collection_name='atoms'):
        self.client = MongoClient(url)
        self.db = self.client[db_name]
        self.atoms = self.db[collection_name]


if __name__ == '__main__':
    from ase.io import iread

    database_name = 'fekad_test'
    collection_name = 'atoms'

    client = MongoClient('mongodb://fekad:qwe123@ds211613.mlab.com:11613/fekad_test')
    db = client[database_name]
    collection = db[collection_name]

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
