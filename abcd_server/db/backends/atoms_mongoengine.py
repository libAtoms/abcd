import logging
import numpy as np
from collections import Counter

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mongoengine import Document, DynamicDocument, EmbeddedDocument, fields, queryset, signals
from mongoengine.context_managers import switch_collection

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


if __name__ == '__main__':
    from mongoengine import connect, register_connection

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
