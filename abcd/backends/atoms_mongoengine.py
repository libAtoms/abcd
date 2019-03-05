import types
import logging
import numpy as np
from os import linesep
from collections import Counter

from ase import Atoms
from ase.io import iread
from ase.calculators.singlepoint import SinglePointCalculator

from mongoengine import Document, DynamicDocument, EmbeddedDocument, fields, queryset, signals, connect
from mongoengine.context_managers import switch_collection

import matplotlib.pyplot as plt

from abcd.backends.base import Database
from abcd.query.parser import QueryParser

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
        else:
            raise NotImplementedError(f'Query string should be string or dictionary and not {type(query)}!')


class MongoQuery(object):

    def __init__(self):
        self.parser = QueryParser()

    def visit(self, syntax_tree):
        op, *args = syntax_tree
        try:
            fun = self.__getattribute__(f'visit_{op.lower()}')
            return fun(*args)
        except:
            pass

    def visit_name(self, fields):
        return {fields: {'$exists': True}}

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
        return {field[1]: value[1]}

    def visit_gt(self, field, value):
        return {field[1]: {'$gt': value[1]}}

    def visit_gte(self, field, value):
        return {field[1]: {'$gte': value[1]}}

    def visit_lt(self, field, value):
        return {field[1]: {'$lt': value[1]}}

    def visit_lte(self, field, value):
        return {field[1]: {'$lte': value[1]}}

    def visit_in(self, field, *values):
        return {field[1]: {'$in': [value[1] for value in values]}}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __call__(self, string):

        ast = self.parser.parse(string)

        if ast is None:
            return {}

        return self.visit(ast)


class DerivedModel(EmbeddedDocument):
    elements = fields.DictField(fields.IntField())
    arrays_keys = fields.ListField(fields.StringField())
    info_keys = fields.ListField(fields.StringField())


class AtomsModel(DynamicDocument):
    meta = {
        'collection': 'atoms',  # default name of the collection
        'queryset_class': AtomsQuerySet
    }

    # arrays = EmbeddedDocumentField(ArraysModel)
    # info = EmbeddedDocumentField(InfoModel)
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
        cell = self['info'].pop('cell', None)
        pbc = self['info'].pop('pbc', None)

        numbers = self['arrays'].pop('numbers', None)
        positions = self['arrays'].pop('positions', None)

        atoms = Atoms(numbers=numbers,
                      cell=cell,
                      pbc=pbc,
                      positions=positions)

        if 'calculator_name' in self['info']:
            calculator_name = self['info'].pop('calculator_name')
            params = self['info'].pop('calculator_parameters', {})
            results = self['results']

            # TODO: Proper initialisation fo Calculators
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            atoms.calc = SinglePointCalculator(atoms, **params, **results)

        atoms.arrays.update(self['arrays'])
        atoms.info.update(self['info'])

        return atoms

    @classmethod
    def from_atoms(cls, atoms: Atoms, extra_info=None, **kwargs):
        arrays = atoms.arrays.copy()
        natoms = len(atoms)

        arrays = {
            'numbers': arrays.pop('numbers').tolist(),
            'positions': arrays.pop('positions').tolist(),

        }
        info = {
            'cell': atoms.cell.tolist(),
            'pbc': atoms.pbc.tolist(),
            'constraints': [],
        }

        for key, value in arrays.items():

            if isinstance(value, np.ndarray):
                arrays[key] = value.tolist()
                continue

            arrays[key] = value

        for key, value in atoms.info.items():

            if isinstance(value, np.ndarray):
                info[key] = value.tolist()
                continue

            info[key] = value

        if atoms.calc is not None:
            info['calculator_name'] = atoms.calc.__class__.__name__
            info['calculator_parameters'] = atoms.calc.todict()

            for key, value in atoms.calc.results.items():

                if isinstance(value, np.ndarray):
                    if value.shape[0] == natoms:
                        arrays[key] = value.tolist()
                    else:
                        info[key] = value.tolist()
                    continue

                info[key] = value

        if extra_info is not None:
            info.update(extra_info)

        return cls(arrays=arrays, info=info)

    @classmethod
    def pre_save_post_validation(cls, sender, document, **kwargs):

        elements = Counter(str(element) for element in document.arrays['numbers'])
        arrays_keys = list(document.arrays.keys())
        info_keys = list(document.info.keys())

        document.derived = DerivedModel(elements=elements, arrays_keys=arrays_keys, info_keys=info_keys)

        logger.debug("Pre Save: %s" % document)

    @classmethod
    def post_save(cls, sender, document, **kwargs):

        logger.debug("Post Save: %s" % document)

        if 'created' in kwargs:
            if kwargs['created']:
                logger.debug("Created")
            else:
                logger.debug("Updated")


signals.pre_save_post_validation.connect(AtomsModel.pre_save_post_validation, sender=AtomsModel)
signals.post_save.connect(AtomsModel.post_save, sender=AtomsModel)


class Databases(Document):
    meta = {'collection': 'databases'}

    name = fields.StringField(required=True)
    description = fields.StringField(max_length=300)
    tags = fields.ListField(fields.StringField(max_length=30))


class MongoDatabase(Database):
    """Wrapper to make database operations easy"""

    def __init__(self, db='abcd', collection='atoms', authentication_source='admin', username=None, password=None,
                 **kwargs):
        super().__init__()

        self.client = connect(db, **kwargs)
        if username:
            self.client[authentication_source].authenticate(name=username, password=password)

        self.collection_name = collection

        self._db = self.client.get_database(db)
        self._collection = self._db[collection]

        self.queryparser = MongoQuery()

    def info(self):
        host, port = self.client.address

        return {
            'host': host,
            'port': port,
            'db': self._db.name,
            'collection': self._collection.name,
            'number of confs': self.count()
        }

    def destroy(self):
        with switch_collection(AtomsModel, self.collection_name) as model:
            model.drop_collection()

    def push(self, atoms, extra_info=None):
        with switch_collection(AtomsModel, self.collection_name) as model:
            if isinstance(atoms, Atoms):
                model.from_atoms(atoms, extra_info).save()

            elif isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):
                # NOTE: insert method is not able to handle generators
                model._get_collection().insert_many(model.from_atoms(at).to_mongo() for at in atoms)
                # model.objects.insert(list(model.from_atoms(at) for at in atoms))

    def upload(self, file):
        data = iread(file)
        self.push(data)

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        raise NotImplementedError

    def get_atoms(self, query=None):
        with switch_collection(AtomsModel, self.collection_name) as model:
            return model.objects.query(query).to_atoms()

    def count(self, query=None):
        with switch_collection(AtomsModel, self.collection_name) as model:
            return model.objects.query(query).count()

    def property(self, name, query=None):
        with switch_collection(AtomsModel, self.collection_name) as model:
            pipeline = [
                {'$match': {'{}'.format(name): {"$exists": True}}},
                {'$project': {'_id': False, 'data': '${}'.format(name)}}
            ]

            return [val['data'] for val in model.objects.query(query).aggregate(*pipeline)]

    def properties(self, query=None):
        with switch_collection(AtomsModel, self.collection_name) as model:
            properties = {}

            pipeline = [
                {'$unwind': '$derived.info_keys'},
                {'$group': {'_id': '$derived.info_keys'}}
            ]
            properties['info'] = [value['_id'] for value in model.objects.query(query).aggregate(*pipeline)]

            pipeline = [
                {'$unwind': '$derived.arrays_keys'},
                {'$group': {'_id': '$derived.arrays_keys'}}
            ]
            properties['arrays'] = [value['_id'] for value in model.objects.query(query).aggregate(*pipeline)]

            return properties

    def count_properties(self, query=None):
        with switch_collection(AtomsModel, self.collection_name) as model:

            properties = {
                'info': {},
                'arrays': {}
            }

            pipeline = [
                {'$unwind': '$derived.info_keys'},
                {'$group': {'_id': '$derived.info_keys', 'count': {'$sum': 1}}}
            ]

            for value in model.objects.query(query).aggregate(*pipeline):
                properties['info'][value['_id']] = {'count': value['count']}

            pipeline = [
                {'$unwind': '$derived.arrays_keys'},
                {'$group': {'_id': '$derived.arrays_keys', 'count': {'$sum': 1}}}
            ]

            for value in model.objects.query(query).aggregate(*pipeline):
                properties['arrays'][value['_id']] = {'count': value['count']}

        return properties

    def __repr__(self):
        host, port = self.client.address

        return f'{self.__class__.__name__}(' \
            f'url={host}:{port}, ' \
            f'db={self._db.name}, ' \
            f'collection={self._collection.name})'

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

    def stats(self, prop, query=None, **kwargs):
        data = self.property(prop, query)
        return {
            prop: {
                'min': min(data),
                'max': max(data),
            }
        }

    def plot_hist(self, prop, query=None, **kwargs):
        data = self.property(prop, query)
        _, _, ax = plt.hist(data, **kwargs)
        return ax


if __name__ == '__main__':
    from ase.io import iread

    logging.basicConfig(level=logging.DEBUG)

    collection_name = 'atoms'
    url = 'mongodb://2ef35d3635e9dc5a922a6a42:ac6ce72e259f5ddcc8dd5178@localhost:27017'

    abcd = MongoDatabase(collection='test', username='2ef35d3635e9dc5a922a6a42', password='ac6ce72e259f5ddcc8dd5178')

    abcd.print_info()

    # for atoms in iread('../../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
    #     # Hack to fix the representation of forces
    #     atoms.calc.results['forces'] = atoms.arrays['force']
    #
    #     print(atoms)

    with switch_collection(AtomsModel, collection_name) as AtomsModel:
        atoms = AtomsModel.objects.first().to_atoms()

    saved = AtomsModel.from_atoms(atoms).save()

    # Update
    Databases(name=collection_name).save()

    print(Databases.objects(name=collection_name).first())

    print(abcd.query('arrays.positions'))

    print(abcd.count_properties())

    query = {
        'info.config_type': 'bcc_bulk_54_high'
    }
    # query ='info.config_type=bcc_bulk_54_high'
    print(abcd.count(query))

# class ArraysModel(EmbeddedDocument):
#     meta = {'strict': False}
#     numbers = ListField(IntField())
#     positions = ListField(ListField(FloatField()))
#
#
# class InfoModel(EmbeddedDocument):
#     meta = {'strict': False}