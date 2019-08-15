import datetime
import getpass
import logging
from collections import Counter

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

logger = logging.getLogger(__name__)


class AbstractModel(dict):

    @classmethod
    def from_atoms(cls, atoms: Atoms):
        """ASE's original implementation"""

        reserved_keys = {'n_atoms', 'cell', 'derived', 'calculator_name', 'calculator_parameters'}
        arrays_keys = set(atoms.arrays.keys())
        info_keys = set(atoms.info.keys())
        results_keys = set(atoms.calc.results.keys()) if atoms.calc else {}

        # all_keys = (reserved_keys, arrays_keys, info_keys, results_keys)
        # if len(set.union(*all_keys)) != sum(map(len, all_keys)):
        #     raise ValueError('All the keys must be unique!')

        n_atoms = len(atoms)

        dct = {
            'n_atoms': n_atoms,
            'cell': atoms.cell.tolist(),
            'pbc': atoms.pbc.tolist(),
        }
        info_keys.update({'n_atoms', 'cell', 'pbc'})

        for key, value in atoms.arrays.items():
            if isinstance(value, np.ndarray):
                dct[key] = value.tolist()
            else:
                dct[key] = value

        for key, value in atoms.info.items():
            if isinstance(value, np.ndarray):
                dct[key] = value.tolist()
            else:
                dct[key] = value

        if atoms.calc is not None:
            dct['calculator_name'] = atoms.calc.__class__.__name__
            dct['calculator_parameters'] = atoms.calc.todict()
            info_keys.update({'calculator_name', 'calculator_parameters'})

            for key, value in atoms.calc.results.items():

                if isinstance(value, np.ndarray):
                    if value.shape[0] == n_atoms:
                        arrays_keys.update(key)
                    else:
                        info_keys.update(key)
                    dct[key] = value.tolist()

        dct['derived'] = {
            'elements': Counter(atoms.get_chemical_symbols()),
            'arrays_keys': list(arrays_keys),
            'info_keys': list(info_keys),
            'results_keys': list(results_keys)
        }

        return cls(**dct)

    def to_atoms(self):
        # TODO: Proper initialisation fo Calculators

        arrays_keys = self['derived']['arrays_keys']
        info_keys = self['derived']['info_keys']
        results_keys = self['derived']['results_keys']

        cell = self.pop('cell', None)
        pbc = self.pop('pbc', None)
        info_keys.remove({'cell', 'pbc'})

        numbers = self.pop('numbers', None)
        arrays_keys.remove('numbers')

        positions = self.pop('positions', None)
        arrays_keys.remove('positions')

        atoms = Atoms(
            numbers=numbers,
            cell=cell,
            pbc=pbc,
            positions=positions)

        if 'calculator_name' in self:
            # calculator_name = self['info'].pop('calculator_name')
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            params = self.pop('calculator_parameters', {})

            atoms.calc = SinglePointCalculator(atoms, **params)
            atoms.calc.results.update((key, self[key]) for key in results_keys)

        atoms.arrays.update((key, self[key]) for key in arrays_keys)
        atoms.arrays.update((key, self[key]) for key in info_keys)

        return atoms

    # # TODO: derived properties when save ()
    # def pre_save_post_validation(cls, sender, document, **kwargs):
    #
    #     document.info['username'] = getpass.getuser()
    #
    #     natoms = len(document.arrays['numbers'])
    #     elements = Counter(str(element) for element in document.arrays['numbers'])
    #
    #     arrays_keys = list(document.arrays.keys())
    #     info_keys = list(document.info.keys())
    #     derived_keys = ['natoms', 'elements', 'username', 'uploaded', 'modified']
    #
    #     cell = document.info.get('cell')
    #     if cell:
    #         derived_keys.append('volume')
    #         virial = document.info.get('virial')
    #         if virial:
    #             derived_keys.append('pressure')
    #
    #     document.derived = dict(
    #         natoms=natoms,
    #         elements=elements,
    #         arrays_keys=arrays_keys,
    #         info_keys=info_keys,
    #         derived_keys=derived_keys,
    #         username=getpass.getuser()
    #     )
    #
    #     cell = document.info.get('cell')
    #     if cell:
    #         volume = abs(np.linalg.det(cell))  # atoms.get_volume()
    #         document.derived['volume'] = volume
    #
    #         virial = document.info.get('virial')
    #         if virial:
    #             # pressure P = -1/3 Tr(stress) = -1/3 Tr(virials/volume)
    #             document.derived['pressure'] = -1 / 3 * np.trace(virial / volume)
    #
    #     if not document.uploaded:
    #         document.uploaded = datetime.datetime.utcnow()
    #
    #     document.modified = datetime.datetime.utcnow()
    #
    #     logger.debug("Pre Save: %s" % document)
    #
    # @classmethod
    # def post_save(cls, sender, document, **kwargs):
    #
    #     logger.debug("Post Save: %s" % document)
    #
    #     if 'created' in kwargs:
    #         if kwargs['created']:
    #             logger.debug("Created")
    #         else:
    #             logger.debug("Updated")
    #
    # @classmethod
    # def pre_bulk_insert(cls, sender, documents, **kwargs):
    #     for document in documents:
    #         cls.pre_save_post_validation(sender, document, **kwargs)


if __name__ == '__main__':
    import io
    from pprint import pprint
    from ase.io import read

    # from ase.io import jsonio

    xyz = io.StringIO("""2
        Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
        Si       0.00000000       0.00000000       0.00000000 
        Si       0.00000000       0.00000000       0.00000000 
        """)

    atoms = read(xyz, format='xyz')
    atoms.set_cell([1, 1, 1])

    # print(atoms)
    # print(atoms.arrays)
    # print(atoms.info)

    # pprint(AbstractModel.from_atoms(atoms))

    # pprint(jsonio.encode(atoms.arrays))
    # pprint(jsonio.encode(atoms.info))
    # pprint(jsonio.encode(atoms.cell))
    #
    pprint(AbstractModel.from_atoms(atoms))

    model = AbstractModel.from_atoms(atoms)
    print(model.to_atoms())
