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
    def __init__(self, url='mongodb://localhost:27017/', db_name=None, collection_name=None):
        self.client = MongoClient(url)
        self.db = self.client[db_name]
        self.atoms = self.db[collection_name]


if __name__ == '__main__':
    from ase.io import iread

    database_name = 'test'
    collection_name = 'atoms'

    client = MongoClient('localhost', 27017)
    db = client[database_name]
    atoms = db[collection_name]

    for at in iread('../../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
        # print(at)
        at.calc.results['forces'] = at.arrays['force']

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
