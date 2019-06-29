import types
import logging
import datetime
import getpass

import numpy as np
from os import linesep
from operator import itemgetter
from collections import Counter

from ase import Atoms
from ase.io import iread
from ase.calculators.singlepoint import SinglePointCalculator

from mongoengine import Document, DynamicDocument, EmbeddedDocument, fields, queryset, signals, connect
import pymongo.errors

from abcd.backends.abstract import Database, URLError, AuthenticationError
from abcd.parsers.queries import QueryParser
from abcd.parsers.arguments import key_val_str_to_dict

logger = logging.getLogger(__name__)


class AtomsQuerySet(queryset.QuerySet):
    """Extension of mongoengine's QuerySet object.
    This class contains helper methods to directly access the Atoms objects
    """

    def from_string(self, query):
        with MongoQuery() as parser:
            query = parser(query)
        return self.from_dict(query)

    def from_dict(self, query):
        return self(__raw__=query)

    def to_atoms(self):
        """Generator which convert all the result of a query to Atoms object.
        """
        return (obj.to_atoms() for obj in self)

    def query(self, query):
        if query is None:
            return self
        elif isinstance(query, dict):
            return self(__raw__=query)
        elif isinstance(query, str):
            return self.from_string(query)
        elif isinstance(query, list):
            return self.from_string(' and '.join(query))
        else:
            raise NotImplementedError('Query string should be string or dictionary and not {}!'.format(type(query)))


class MongoQuery(object):

    def __init__(self):
        self.parser = QueryParser()

    def visit(self, syntax_tree):
        op, *args = syntax_tree
        try:
            fun = self.__getattribute__('visit_' + op.lower())
            return fun(*args)
        except:
            pass

    def visit_name(self, fields):
        # return {fields: {'$exists': True}}
        return {
            '$or': [
                {'info.' + fields: {'$exists': True}},
                {'arrays.' + fields: {'$exists': True}},
                {'derived.' + fields: {'$exists': True}},
            ]
        }

    def visit_and(self, *args):
        return {'$and': [self.visit(arg) for arg in args]}
        # TODO recursively combining all the and statements
        # out = {}
        # for arg in args:
        #     a = self.visit(arg)
        #
        #     out.update(**a)
        # return out

    def visit_or(self, *args):
        return {'$or': [self.visit(arg) for arg in args]}

    def visit_eq(self, field, value):
        # return {field[1]: value[1]}
        return {
            '$or': [
                {'info.' + field[1]: value[1]},
                {'arrays.' + field[1]: value[1]},
                {'derived.' + field[1]: value[1]},
            ]
        }

    def visit_re(self, field, value):
        # return {field[1]: {'$regex': value[1]}}
        return {
            '$or': [
                {'info.' + field[1]: {'$regex': value[1]}},
                {'arrays.' + field[1]: {'$regex': value[1]}},
                {'derived.' + field[1]: {'$regex': value[1]}},
            ]
        }

    def visit_gt(self, field, value):
        # return {field[1]: {'$gt': value[1]}}
        return {
            '$or': [
                {'info.' + field[1]: {'$gt': value[1]}},
                {'arrays.' + field[1]: {'$gt': value[1]}},
                {'derived.' + field[1]: {'$gt': value[1]}},
            ]
        }

    def visit_gte(self, field, value):
        # return {field[1]: {'$gte': value[1]}}
        return {
            '$or': [
                {'info.' + field[1]: {'$gte': value[1]}},
                {'arrays.' + field[1]: {'$gte': value[1]}},
                {'derived.' + field[1]: {'$gte': value[1]}},
            ]
        }

    def visit_lt(self, field, value):
        # return {field[1]: {'$lt': value[1]}}
        return {
            '$or': [
                {'info.' + field[1]: {'$lt': value[1]}},
                {'arrays.' + field[1]: {'$lt': value[1]}},
                {'derived.' + field[1]: {'$lt': value[1]}},
            ]
        }

    def visit_lte(self, field, value):
        # return {field[1]: {'$lte': value[1]}}
        return {
            '$or': [
                {'info.' + field[1]: {'$lte': value[1]}},
                {'arrays.' + field[1]: {'$lte': value[1]}},
                {'derived.' + field[1]: {'$lte': value[1]}},
            ]
        }

    def visit_in(self, field, *values):
        # return {field[1]: {'$in': [value[1] for value in values]}}
        return {
            '$or': [
                {'info.' + field[1]: {'$in': [value[1] for value in values]}},
                {'arrays.' + field[1]: {'$in': [value[1] for value in values]}},
                {'derived.' + field[1]: {'$in': [value[1] for value in values]}},
            ]
        }

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __call__(self, string):
        ast = self.parser.parse(string)
        return self.visit(ast) if ast is not None else {}


class DerivedModel(EmbeddedDocument):
    username = fields.StringField()
    natoms = fields.IntField()
    elements = fields.DictField(fields.IntField())
    volume = fields.FloatField()
    pressure = fields.FloatField()

    arrays_keys = fields.ListField(fields.StringField())
    info_keys = fields.ListField(fields.StringField())
    derived_keys = fields.ListField(fields.StringField())


class AtomsModel(DynamicDocument):
    meta = {
        'collection': 'atoms',  # default name of the collection
        'queryset_class': AtomsQuerySet
    }

    # arrays = EmbeddedDocumentField(ArraysModel)
    # info = EmbeddedDocumentField(InfoModel)

    author = fields.StringField(max_length=120, default='anonymous')
    uploaded = fields.DateTimeField()
    modified = fields.DateTimeField()

    arrays = fields.DictField()
    info = fields.DictField()
    results = fields.DictField()
    derived = fields.EmbeddedDocumentField(DerivedModel)

    def to_atoms(self, *args, **kwargs):
        """Convert the query to ASE Atoms object.

        Args:
            *args:
            **kwargs:

        Returns:

        """
        # TODO: properly converting lists to numpy arrays

        cell = self['info'].pop('cell', None)
        pbc = self['info'].pop('pbc', None)

        numbers = self['arrays'].pop('numbers', None)
        positions = self['arrays'].pop('positions', None)

        atoms = Atoms(numbers=np.array(numbers),
                      cell=np.array(cell),
                      pbc=np.array(pbc),
                      positions=np.array(positions))

        if 'calculator_name' in self['info']:
            calculator_name = self['info'].pop('calculator_name')
            params = self['info'].pop('calculator_parameters', {})
            results = self['results']

            for k, v in results.items():
                results[k] = np.array(v)

            # TODO: Proper initialisation fo Calculators
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            atoms.calc = SinglePointCalculator(atoms, **params, **results)

        for k, v in self['arrays'].items():
            self['arrays'][k] = np.array(v)

        atoms.arrays.update(self['arrays'])
        atoms.info.update(self['info'])

        return atoms

    @classmethod
    def from_atoms(cls, atoms: Atoms, extra_info=None, **kwargs):

        # if isinstance(atoms, list):
        #     return (cls.from_atoms(value, extra_info=None, **kwargs) for value in atoms)
        arrays = atoms.arrays.copy()
        natoms = len(atoms)

        db_arrays = {
            'numbers': arrays.pop('numbers').tolist(),
            'positions': arrays.pop('positions').tolist(),
        }
        db_info = {
            'cell': atoms.cell.tolist(),
            'pbc': atoms.pbc.tolist(),
            'constraints': [],
            'formula': atoms.get_chemical_formula()
        }

        for key, value in arrays.items():

            if isinstance(value, np.ndarray):
                db_arrays[key] = value.tolist()
                continue

            db_arrays[key] = value

        for key, value in atoms.info.items():

            if isinstance(value, np.ndarray):
                db_info[key] = value.tolist()
                continue
            elif isinstance(value, np.int64):
                db_info[key] = int(value)
                continue
            elif isinstance(value, np.float64):
                db_info[key] = float(value)
                continue

            db_info[key] = value

        if atoms.calc is not None:
            db_info['calculator_name'] = atoms.calc.__class__.__name__
            db_info['calculator_parameters'] = atoms.calc.todict()

            for key, value in atoms.calc.results.items():

                if isinstance(value, np.ndarray):
                    if value.shape[0] == natoms:
                        db_arrays[key] = value.tolist()
                    else:
                        db_info[key] = value.tolist()
                    continue

                db_info[key] = value

        if extra_info is not None:
            db_info.update(extra_info)

        return cls(arrays=db_arrays, info=db_info)

    @classmethod
    def pre_save_post_validation(cls, sender, document, **kwargs):

        document.info['username'] = getpass.getuser()

        natoms = len(document.arrays['numbers'])
        elements = Counter(str(element) for element in document.arrays['numbers'])

        arrays_keys = list(document.arrays.keys())
        info_keys = list(document.info.keys())
        derived_keys = ['natoms', 'elements', 'username', 'uploaded', 'modified']

        cell = document.info.get('cell')
        if cell:
            derived_keys.append('volume')
            virial = document.info.get('virial')
            if virial:
                derived_keys.append('pressure')

        document.derived = DerivedModel(
            natoms=natoms,
            elements=elements,
            arrays_keys=arrays_keys,
            info_keys=info_keys,
            derived_keys=derived_keys,
            username=getpass.getuser()
        )

        cell = document.info.get('cell')
        if cell:
            volume = abs(np.linalg.det(cell))  # atoms.get_volume()
            document.derived['volume'] = volume

            virial = document.info.get('virial')
            if virial:
                # pressure P = -1/3 Tr(stress) = -1/3 Tr(virials/volume)
                document.derived['pressure'] = -1 / 3 * np.trace(virial / volume)

        if not document.uploaded:
            document.uploaded = datetime.datetime.utcnow()

        document.modified = datetime.datetime.utcnow()

        logger.debug("Pre Save: %s" % document)

    @classmethod
    def post_save(cls, sender, document, **kwargs):

        logger.debug("Post Save: %s" % document)

        if 'created' in kwargs:
            if kwargs['created']:
                logger.debug("Created")
            else:
                logger.debug("Updated")

    @classmethod
    def pre_bulk_insert(cls, sender, documents, **kwargs):
        for document in documents:
            cls.pre_save_post_validation(sender, document, **kwargs)


signals.pre_save_post_validation.connect(AtomsModel.pre_save_post_validation, sender=AtomsModel)
signals.pre_bulk_insert.connect(AtomsModel.pre_bulk_insert, sender=AtomsModel)

signals.post_save.connect(AtomsModel.post_save, sender=AtomsModel)


class Databases(Document):
    meta = {'collection': 'databases'}

    name = fields.StringField(required=True)
    description = fields.StringField(max_length=300)
    tags = fields.ListField(fields.StringField(max_length=30))


class MongoDatabase(Database):
    """Wrapper to make database operations easy"""

    def __init__(self, db='abcd', collection='atoms',
                 username=None, password=None, authentication_source='admin',
                 **kwargs):
        super().__init__()

        try:
            self.client = connect(db, **kwargs)

            if username:
                self.client[authentication_source].authenticate(name=username, password=password)

        except pymongo.errors.ServerSelectionTimeoutError:
            raise URLError() from None
        except pymongo.errors.OperationFailure:
            raise AuthenticationError() from None

        self.db_name = db
        self.collection_name = collection

        self._db = self.client.get_database(db)
        self._collection = self._db[collection]

        self.queryparser = MongoQuery()

    def info(self):
        host, port = self.client.address

        return {
            'host': host,
            'port': port,
            'db': self.db_name,
            'collection': self.collection_name,
            'number of confs': self.count(),
            'type': 'mongodb'
        }

    def delete(self, query=None):
        return AtomsModel.objects.query(query).delete()

    def destroy(self):
        AtomsModel.drop_collection()

    def push(self, atoms, extra_info=None):

        if extra_info and isinstance(extra_info, str):
            extra_info = key_val_str_to_dict(extra_info)

        if isinstance(atoms, Atoms):
            AtomsModel.from_atoms(atoms, extra_info).save()

        elif isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):
            # NOTE: insert method is not able to handle generators

            # model._get_collection().insert_many(model.from_atoms(at).to_mongo() for at in atoms)
            AtomsModel.objects.insert(list(AtomsModel.from_atoms(at, extra_info) for at in atoms))

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

    def get_atoms(self, query=None):
        return AtomsModel.objects.query(query).to_atoms()

    def count(self, query=None):
        return AtomsModel.objects.query(query).count()

    def property(self, name, query=None):
        pipeline = [
            {'$match': {'{}'.format(name): {"$exists": True}}},
            {'$project': {'_id': False, 'data': '${}'.format(name)}}
        ]

        return [val['data'] for val in AtomsModel.objects.query(query).aggregate(*pipeline)]

    def properties(self, query=None):
        properties = {}

        pipeline = [
            {'$unwind': '$derived.info_keys'},
            {'$group': {'_id': '$derived.info_keys'}}
        ]
        properties['info'] = [value['_id'] for value in AtomsModel.objects.query(query).aggregate(*pipeline)]

        pipeline = [
            {'$unwind': '$derived.arrays_keys'},
            {'$group': {'_id': '$derived.arrays_keys'}}
        ]
        properties['arrays'] = [value['_id'] for value in AtomsModel.objects.query(query).aggregate(*pipeline)]

        return properties

    def count_properties(self, query=None):

        properties = {
            'info': {},
            'arrays': {},
            'derived': {}
        }

        pipeline = [
            {'$unwind': '$derived.info_keys'},
            {'$group': {'_id': '$derived.info_keys', 'count': {'$sum': 1}}}
        ]

        for value in AtomsModel.objects.query(query).aggregate(*pipeline):
            properties['info'][value['_id']] = {'count': value['count']}

        pipeline = [
            {'$unwind': '$derived.arrays_keys'},
            {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}
        ]

        for value in AtomsModel.objects.query(query).aggregate(*pipeline):
            properties['arrays'][value['_id']] = {'count': value['count']}

        pipeline = [
            {'$unwind': '$derived.derived_keys'},
            {'$group': {'_id': '$derived.derived_keys', 'count': {'$sum': 1}}}
        ]

        for value in AtomsModel.objects.query(query).aggregate(*pipeline):
            properties['derived'][value['_id']] = {'count': value['count']}

        return properties

    # def rename_properties(self, name, new_name, query=None, overwrite=False):
    #     print(name, new_name, query, overwrite)

    def __repr__(self):
        host, port = self.client.address

        return '{}('.format(self.__class__.__name__) + \
               'url={}:{}, '.format(host, port) + \
               'db={}, '.format(self.db_name) + \
               'collection={})'.format(self.collection_name)

    def _repr_html_(self):
        """Jupyter notebook representation"""
        return '<b>ABCD MongoDB database</b>'

    # noinspection PyStringFormat
    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join(
            ['{:=^50}'.format(' ABCD MongoDB '),
             '{:>10}: {}'.format('type', 'mongodb'),
             linesep.join('{:>10}: {}'.format(k, v) for k, v in self.info().items())])

        print(out)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def stats(self, prop, query=None, **kwargs):
        data = self.property(prop, query)
        return {
            prop: {
                'min': min(data),
                'max': max(data),
            }
        }

    def hist(self, name, query=None, **kwargs):

        data = self.property(name, query)

        if not data:
            return None

        elif data and isinstance(data, list):
            if isinstance(data[0], float):
                bins = kwargs.get('bins', 10)
                return self._hist_float(name, data, bins)

            elif isinstance(data[0], int):
                bins = kwargs.get('bins', 10)
                return self._hist_int(name, data, bins)

            elif isinstance(data[0], str):
                return self._hist_str(name, data, **kwargs)

            else:
                print('{}: Histogram for list of {} types are not supported!'.format(name, type(data)))
                logger.info('{}: Histogram for list of {} types are not supported!'.format(name, type(data)))

        else:
            logger.info('{}: Histogram for {} types are not supported!'.format(name, type(data)))
            return None

    @staticmethod
    def _hist_float(name, data, bins=10):

        data = np.array(data)
        hist, bin_edges = np.histogram(data, bins=bins)

        return {
            'type': 'hist_float',
            'name': name,
            'bins': bins,
            'edges': bin_edges,
            'counts': hist,
            'min': data.min(),
            'max': data.max(),
            'median': data.mean(),
            'std': data.std(),
            'var': data.var()
        }

    @staticmethod
    def _hist_int(name, data, bins=10):

        data = np.array(data)
        delta = max(data) - min(data) + 1

        if bins > delta:
            bins = delta

        hist, bin_edges = np.histogram(data, bins=bins)

        return {
            'type': 'hist_int',
            'name': name,
            'bins': bins,
            'edges': bin_edges,
            'counts': hist,
            'min': data.min(),
            'max': data.max(),
            'median': data.mean(),
            'std': data.std(),
            'var': data.var()
        }

    @staticmethod
    def _hist_str(name, data, bins=10, truncate=20):

        if truncate:
            # data = (item[:truncate] for item in data)
            data = (item[:truncate] + '...' if len(item) > truncate else item for item in data)

        data = Counter(data)

        if bins:
            labels, counts = zip(*sorted(data.items(), key=itemgetter(1, 0), reverse=True))
        else:
            labels, counts = zip(*data.items())

        return {
            'type': 'hist_str',
            'name': name,
            'total': sum(data.values()),
            'unique': len(data.keys()),
            'labels': labels[:bins],
            'counts': counts[:bins]
        }

    def exec(self, code, query=None):
        for item in AtomsModel.objects.query(query):
            exec(code)


# def debug_issue19():
#     from pathlib import Path
#     from ase.io import read
#
#     file = Path('../../gp_iter6_sparse9k.xml.xyz')
#     traj = read(file.as_posix(), index='170:172')
#
#     abcd = MongoDatabase()
#     abcd.push(traj[1])
#     abcd.push(traj)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    # url = 'mongodb://mongoadmin:secret@localhost:27017/abcd'
    abcd = MongoDatabase(db='abcd', username='mongoadmin', password='secret')
    # abcd.exec('print(item.info["energy"])', query='natoms=123')

    #
    # def exp(at):
    #     at.info['pressure'] = -1/3 * np.trace(at.info['virial']/at.get_volume())
    #
    # {
    #     'pressure': lambda at: -1/3 * np.trace(at.info['virial']/at.get_volume())
    # }

    # debug_issue19()
#
#     collection_name = 'atoms'
#     url = 'mongodb://2ef35d3635e9dc5a922a6a42:ac6ce72e259f5ddcc8dd5178@localhost:27017'
#
#     abcd = MongoDatabase(collection='atoms', username='2ef35d3635e9dc5a922a6a42', password='ac6ce72e259f5ddcc8dd5178')
#
#     abcd.print_info()
#
#     abcd.push(iread('../../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(10)))
#     # for atoms in iread('../../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
#     #     # Hack to fix the representation of forces
#     #     atoms.calc.results['forces'] = atoms.arrays['force']
#     #
#     #     print(atoms)
#
#     # with switch_collection(AtomsModel, collection_name) as AtomsModel:
#     atoms = AtomsModel.objects.first().to_atoms()
#
#     saved = AtomsModel.from_atoms(atoms).save()
#
#     # Update
#     Databases(name=collection_name).save()
#
#     print(Databases.objects(name=collection_name).first())
#
#     print(abcd.query('arrays.positions'))
#
#     print(abcd.count_properties())
#
#     query = {
#         'info.config_type': 'bcc_bulk_54_high'
#     }
#     # query ='info.config_type=bcc_bulk_54_high'
#     print(abcd.count(query))
#
#     a = list(abcd.get_atoms())
#     print(type(a[0].arrays['forces']))
# # class ArraysModel(EmbeddedDocument):
# #     meta = {'strict': False}
# #     numbers = ListField(IntField())
# #     positions = ListField(ListField(FloatField()))
# #
# #
# # class InfoModel(EmbeddedDocument):
# #     meta = {'strict': False}
