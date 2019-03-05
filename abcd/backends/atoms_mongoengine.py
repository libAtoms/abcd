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

from abcd.backends.base import Database

logger = logging.getLogger(__name__)


class AtomsQuerySet(queryset.QuerySet):
    """Extension of mongoengine's QuerySet object.
    This class contains helper methods to directly access the Atoms objects
    """

    def to_atoms(self):
        """Generator which convert all the result of a query to Atoms object.
        """
        # return [obj.to_atoms() for obj in self]

        for obj in self:
            yield obj.to_atoms()


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

    def __init__(self, db='abcd', collection='atoms',  **kwargs):
        super().__init__()

        self.client = connect(db, **kwargs)
        # source='admin' is the default...
        self.client[kwargs['authentication_source']].authenticate(
            name=kwargs['username'], password=kwargs['password']
        )

        self.collection_name = collection

        self.db = self.client.get_database(db)
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

    def query(self, query_string):
        raise NotImplementedError

    def search(self, query_string: str):
        raise NotImplementedError

    def get_atoms(self, dbfilter):
        for dct in self.db.atoms.find(dbfilter):
            yield AtomsModel(dct).to_atoms()

    def count(self, dbfilter=None):
        with switch_collection(AtomsModel, self.collection_name) as model:
            if dbfilter is None:
                return model.objects.count()

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
    logging.basicConfig(level=logging.DEBUG)

    collection_name = 'atoms'

    # register_connection('default', db='abcd', host='localhost', port=27017)
    con = connect('abcd', host='localhost', port=27017)
    db = con.get_database('abcd')
    print(db)

    with switch_collection(AtomsModel, collection_name) as AtomsModel:
        atoms = AtomsModel.objects.first().to_atoms()

    saved = AtomsModel.from_atoms(atoms).save()

    # Update
    Databases(name=collection_name).save()

    print(Databases.objects(name=collection_name).first())

# class ArraysModel(EmbeddedDocument):
#     meta = {'strict': False}
#     numbers = ListField(IntField())
#     positions = ListField(ListField(FloatField()))
#
#
# class InfoModel(EmbeddedDocument):
#     meta = {'strict': False}
