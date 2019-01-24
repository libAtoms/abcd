from abcd_server.encoders.base import BaseEncoder
import numpy as np


# from abc import ABCMeta, abstractmethod
#
#
# class BaseEncoder(object, metaclass=ABCMeta):
#     """Abstract class for the visitor pattern"""
#     default_properties = []
#
#     # @abstractmethod
#     def __init__(self):
#         pass
#
#     def __enter__(self):
#         """support with statement and error handling in python"""
#         return self
#
#     def __exit__(self, exc_type, exc_val, exc_tb):
#         pass
#
#     def encode(self, atoms):
#         """main function"""
#         return self.visit_atoms(atoms)
#
#     @abstractmethod
#     def visit_atoms(self, atoms):
#         pass
#
#     @abstractmethod
#     def visit_numbers(self, atoms):
#         pass
#
#     @abstractmethod
#     def visit_cell(self, atoms):
#         pass
#
#     @abstractmethod
#     def visit_pbc(self, atoms):
#         pass
#
#     @abstractmethod
#     def visit_positions(self, atoms):
#         pass
#
#     @abstractmethod
#     def visit_forces(self, atoms):
#         pass
#
#     @abstractmethod
#     def visit_energy(self, atoms):
#         pass


class DictEncoder(BaseEncoder):
    # default_properties = ['numbers', 'positions']

    def __init__(self):
        super().__init__()

    def visit_atoms(self, atoms):
        """walk through on the structure of an Atoms object """
        data = {
            'numbers': self.visit_numbers(atoms),
            'cell': self.visit_cell(atoms),
            'pbc': self.visit_pbc(atoms),
            'positions': self.visit_positions(atoms),
            'forces': self.visit_forces(atoms),
            'energy': self.visit_energy(atoms),
            'metadata': self.visit_metadata(atoms)
        }
        return data

    def visit_metadata(self, atoms):
        # metadata = {}
        # for key, val in atoms.arrays.items():
        #     if key not in self.default_properties:
        #         metadata[key] = self.convert(val)
        return {key: self.convert(val) for key, val in atoms.arrays.items() if key not in self.default_properties}

    def visit_numbers(self, atoms):
        return super().visit_numbers(atoms).tolist()

    def visit_cell(self, atoms):
        return super().visit_cell(atoms).tolist()

    def visit_pbc(self, atoms):
        return super().visit_pbc(atoms).tolist()

    def visit_positions(self, atoms):
        return super().visit_positions(atoms).tolist()

    def visit_forces(self, atoms):
        return super().visit_forces(atoms).tolist()

    def convert(self, data):
        if isinstance(data, np.ndarray):
            return data.tolist()
        if isinstance(data, dict):
            return {key: self.convert(val) for key, val in data.items()}
        return data


class OrderedDictEncoder(BaseEncoder):
    default_properties = ['numbers', 'positions']

    def __init__(self):
        super().__init__()

    def visit_atoms(self, atoms):
        """walk through on the structure of an Atoms object """
        data = OrderedDict([
            ('unique_id', '{}'.format(randint(16 ** 31, 16 ** 32 - 1))),
            ('numbers', self.visit_numbers(atoms)),
            ('positions', self.visit_positions(atoms)),
            ('pbc', self.visit_pbc(atoms)),
            ('cell', self.visit_cell(atoms)),
            # ('initial_magmoms', self.visit_initial_magmoms(atoms)),
            # ('initial_charges', self.visit_initial_charges(atoms)),
            # ('masses', self.visit_masses(atoms)),
            # ('tags', self.visit_tags(atoms)),
            # ('momenta', self.visit_momenta(atoms)),
            # ('constraints', self.visit_constraints(atoms)),
            # ('calculator', self.visit_calculator(atoms)),
        ])
        return data

    def visit_metadata(self, atoms):
        # metadata = {}
        # for key, val in atoms.arrays.items():
        #     if key not in self.default_properties:
        #         metadata[key] = self.convert(val)
        return {key: self.convert(val) for key, val in atoms.arrays.items() if key not in self.default_properties}

    def visit_numbers(self, atoms):
        return super().visit_numbers(atoms).tolist()

    def visit_cell(self, atoms):
        return super().visit_cell(atoms).tolist()

    def visit_pbc(self, atoms):
        return super().visit_pbc(atoms).tolist()

    def visit_positions(self, atoms):
        return super().visit_positions(atoms).tolist()

    def visit_forces(self, atoms):
        return super().visit_forces(atoms).tolist()

    def visit_initial_magmoms(self, atoms):
        return atoms.get_initial_magmoms()

    def visit_initial_charges(self, atoms):
        return atoms.get_initial_charges()

    def visit_masses(self, atoms):
        return atoms.get_masses()

    def visit_tags(self, atoms):
        return atoms.get_tags()

    def visit_momenta(self, atoms):
        return atoms.get_momenta()

    def visit_constraints(self, atoms):
        return atoms.get_constraints()

    def visit_calculator(self, atoms):
        return atoms.get_calculator()

    def convert(self, data):
        if isinstance(data, np.ndarray):
            return data.tolist()
        if isinstance(data, dict):
            return {key: self.convert(val) for key, val in data.items()}
        return data


# from ase


import json
import datetime
from random import randint
from ase import Atoms
from ase.utils import basestring, formula_metal
from ase.data import chemical_symbols, atomic_masses

from ase.calculators.calculator import get_calculator
from ase.calculators.singlepoint import SinglePointCalculator

import numpy as np

all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy']

all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']

# Recognized names of calculators sorted alphabetically:
names = ['abinit', 'aims', 'amber', 'asap', 'castep', 'cp2k', 'crystal',
         'demon', 'dftb', 'dmol', 'eam', 'elk', 'emt', 'espresso',
         'exciting', 'fleur', 'gaussian', 'gpaw', 'gromacs', 'gulp',
         'hotbit', 'jacapo', 'lammpsrun',
         'lammpslib', 'lj', 'mopac', 'morse', 'nwchem', 'octopus', 'onetep',
         'openmx', 'siesta', 'tip3p', 'turbomole', 'vasp']

special = {'cp2k': 'CP2K',
           'dmol': 'DMol3',
           'eam': 'EAM',
           'elk': 'ELK',
           'emt': 'EMT',
           'crystal': 'CRYSTAL',
           'fleur': 'FLEUR',
           'gulp': 'GULP',
           'lammpsrun': 'LAMMPS',
           'lammpslib': 'LAMMPSlib',
           'lj': 'LennardJones',
           'mopac': 'MOPAC',
           'morse': 'MorsePotential',
           'nwchem': 'NWChem',
           'openmx': 'OpenMX',
           'tip3p': 'TIP3P'}


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


# def get_calculator(name):
#     """Return calculator class."""
#     if name == 'asap':
#         from asap3 import EMT as Calculator
#     elif name == 'gpaw':
#         from gpaw import GPAW as Calculator
#     elif name == 'hotbit':
#         from hotbit import Calculator
#     elif name == 'vasp2':
#         from ase.calculators.vasp import Vasp2 as Calculator
#     else:
#         classname = special.get(name, name.title())
#         module = __import__('ase.calculators.' + name, {}, None, [classname])
#         Calculator = getattr(module, classname)
#     return Calculator

# TODO: OrderedDict!
# TODO: if has that attribute the get function is completely unnecessary
# TODO: cell and pbc is always presents
# TODO: has in a redundant function call
# TODO: calculator_parameters should be an OrderedDict
# TODO: unique id should be a hash value

from collections import OrderedDict


def to_dict(atoms: Atoms):
    dct = OrderedDict([
        ('numbers', atoms.numbers),
        ('positions', atoms.positions),
        ('unique_id', '{}'.format(randint(16 ** 31, 16 ** 32 - 1)))
    ])

    if atoms.cell.any():
        dct['pbc'] = atoms.pbc
        dct['cell'] = atoms.cell
    if atoms.has('initial_magmoms'):
        dct['initial_magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.has('initial_charges'):
        dct['initial_charges'] = atoms.get_initial_charges()
    if atoms.has('masses'):
        dct['masses'] = atoms.get_masses()
    if atoms.has('tags'):
        dct['tags'] = atoms.get_tags()
    if atoms.has('momenta'):
        dct['momenta'] = atoms.get_momenta()
    if atoms.constraints:
        dct['constraints'] = [c.todict() for c in atoms.constraints]
    if atoms.calc is not None:
        dct['calculator'] = atoms.calc.name.lower()
        dct['calculator_parameters'] = atoms.calc.todict()
        if len(atoms.calc.check_state(atoms)) == 0:
            for prop in all_properties:
                try:
                    x = atoms.calc.get_property(prop, atoms, False)
                except PropertyNotImplementedError:
                    pass
                else:
                    if x is not None:
                        dct[prop] = x
    return dct


def atoms2dict(atoms):
    """ASE's original implementation"""
    dct = {
        'numbers': atoms.numbers,
        'positions': atoms.positions,
        'unique_id': '{}'.format(randint(16 ** 31, 16 ** 32 - 1))
    }
    if atoms.cell.any():
        dct['pbc'] = atoms.pbc
        dct['cell'] = atoms.cell
    if atoms.has('initial_magmoms'):
        dct['initial_magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.has('initial_charges'):
        dct['initial_charges'] = atoms.get_initial_charges()
    if atoms.has('masses'):
        dct['masses'] = atoms.get_masses()
    if atoms.has('tags'):
        dct['tags'] = atoms.get_tags()
    if atoms.has('momenta'):
        dct['momenta'] = atoms.get_momenta()
    if atoms.constraints:
        dct['constraints'] = [c.todict() for c in atoms.constraints]
    if atoms.calc is not None:
        dct['calculator'] = atoms.calc.name.lower()
        dct['calculator_parameters'] = atoms.calc.todict()
        if len(atoms.calc.check_state(atoms)) == 0:
            for prop in all_properties:
                try:
                    x = atoms.calc.get_property(prop, atoms, False)
                except PropertyNotImplementedError:
                    pass
                else:
                    if x is not None:
                        dct[prop] = x
    return dct


def numpyfy(obj):
    if isinstance(obj, dict):
        if '__complex_ndarray__' in obj:
            r, i = (np.array(x) for x in obj['__complex_ndarray__'])
            return r + i * 1j
        return dict((intkey(key), numpyfy(value))
                    for key, value in obj.items())
    if isinstance(obj, list) and len(obj) > 0:
        try:
            a = np.array(obj)
        except ValueError:
            pass
        else:
            if a.dtype in [bool, int, float]:
                return a
        obj = [numpyfy(value) for value in obj]
    return obj


def intkey(key):
    try:
        return int(key)
    except ValueError:
        return key


def object_hook(dct):
    if '__datetime__' in dct:
        return datetime.datetime.strptime(dct['__datetime__'],
                                          '%Y-%m-%dT%H:%M:%S.%f')
    if '__complex_ndarray__' in dct:
        r, i = (np.array(x) for x in dct['__complex_ndarray__'])
        return r + i * 1j
    return dct


mydecode = json.JSONDecoder(object_hook=object_hook).decode


def decode(txt):
    return numpyfy(mydecode(txt))


my__all__ = ['FixCartesian', 'FixBondLength', 'FixedMode', 'FixConstraintSingle',
             'FixAtoms', 'UnitCellFilter', 'ExpCellFilter', 'FixScaled', 'StrainFilter',
             'FixCom', 'FixedPlane', 'Filter', 'FixConstraint', 'FixedLine',
             'FixBondLengths', 'FixInternals', 'Hookean', 'ExternalForce']


def dict2constraint(dct):
    if dct['name'] not in my__all__:
        raise ValueError
    return globals()[dct['name']](**dct['kwargs'])


class AtomsRow:
    def __init__(self, dct):
        dct = dct.copy()
        if 'calculator_parameters' in dct:
            # Earlier version of ASE would encode the calculator
            # parameter dict again and again and again ...
            while isinstance(dct['calculator_parameters'], basestring):
                dct['calculator_parameters'] = decode(
                    dct['calculator_parameters'])

        assert 'numbers' in dct
        self._constraints = dct.pop('constraints', [])
        self._constrained_forces = None
        self._data = dct.pop('data', {})
        kvp = dct.pop('key_value_pairs', {})
        self._keys = list(kvp.keys())
        self.__dict__.update(kvp)
        self.__dict__.update(dct)
        if 'cell' not in dct:
            self.cell = np.zeros((3, 3))
            self.pbc = np.zeros(3, bool)

    def __contains__(self, key):
        return key in self.__dict__

    def __iter__(self):
        return (key for key in self.__dict__ if key[0] != '_')

    def get(self, key, default=None):
        """Return value of key if present or default if not."""
        return getattr(self, key, default)

    @property
    def key_value_pairs(self):
        """Return dict of key-value pairs."""
        return dict((key, self.get(key)) for key in self._keys)

    def count_atoms(self):
        """Count atoms.

        Return dict mapping chemical symbol strings to number of atoms.
        """
        count = {}
        for symbol in self.symbols:
            count[symbol] = count.get(symbol, 0) + 1
        return count

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __str__(self):
        return '<AtomsRow: formula={0}, keys={1}>'.format(
            self.formula, ','.join(self._keys))

    @property
    def constraints(self):
        """List of constraints."""
        if not isinstance(self._constraints, list):
            # Lazy decoding:
            cs = decode(self._constraints)
            self._constraints = []
            for c in cs:
                # Convert to new format:
                name = c.pop('__name__', None)
                if name:
                    c = {'name': name, 'kwargs': c}
                if c['name'].startswith('ase'):
                    c['name'] = c['name'].rsplit('.', 1)[1]
                self._constraints.append(c)
        return [dict2constraint(d) for d in self._constraints]

    @property
    def data(self):
        """Data dict."""
        if not isinstance(self._data, dict):
            self._data = decode(self._data)  # lazy decoding
        return FancyDict(self._data)

    @property
    def natoms(self):
        """Number of atoms."""
        return len(self.numbers)

    @property
    def formula(self):
        """Chemical formula string."""
        return formula_metal(self.numbers)

    @property
    def symbols(self):
        """List of chemical symbols."""
        return [chemical_symbols[Z] for Z in self.numbers]

    @property
    def fmax(self):
        """Maximum atomic force."""
        forces = self.constrained_forces
        return (forces ** 2).sum(1).max() ** 0.5

    @property
    def constrained_forces(self):
        """Forces after applying constraints."""
        if self._constrained_forces is not None:
            return self._constrained_forces
        forces = self.forces
        constraints = self.constraints
        if constraints:
            forces = forces.copy()
            atoms = self.toatoms()
            for constraint in constraints:
                constraint.adjust_forces(atoms, forces)

        self._constrained_forces = forces
        return forces

    @property
    def smax(self):
        """Maximum stress tensor component."""
        return (self.stress ** 2).max() ** 0.5

    @property
    def mass(self):
        """Total mass."""
        if 'masses' in self:
            return self.masses.sum()
        return atomic_masses[self.numbers].sum()

    @property
    def volume(self):
        """Volume of unit cell."""
        if self.cell is None:
            return None
        vol = abs(np.linalg.det(self.cell))
        if vol == 0.0:
            raise AttributeError
        return vol

    @property
    def charge(self):
        """Total charge."""
        charges = self.get('inital_charges')
        if charges is None:
            return 0.0
        return charges.sum()

    def toatoms(self, attach_calculator=False,
                add_additional_information=False):
        """Create Atoms object."""
        atoms = Atoms(self.numbers,
                      self.positions,
                      cell=self.cell,
                      pbc=self.pbc,
                      magmoms=self.get('initial_magmoms'),
                      charges=self.get('initial_charges'),
                      tags=self.get('tags'),
                      masses=self.get('masses'),
                      momenta=self.get('momenta'),
                      constraint=self.constraints)

        if attach_calculator:
            params = self.get('calculator_parameters', {})
            atoms.calc = get_calculator(self.calculator)(**params)
        else:
            results = {}
            for prop in all_properties:
                if prop in self:
                    results[prop] = self[prop]
            if results:
                atoms.calc = SinglePointCalculator(atoms, **results)
                atoms.calc.name = self.get('calculator', 'unknown')

        if add_additional_information:
            atoms.info = {}
            atoms.info['unique_id'] = self.unique_id
            if self._keys:
                atoms.info['key_value_pairs'] = self.key_value_pairs
            data = self.get('data')
            if data:
                atoms.info['data'] = data

        return atoms


from ase.io.extxyz import *

if __name__ == '__main__':
    from pathlib import Path
    from ase.io import iread

    from ase.db.row import AtomsRow, atoms2dict

    direcotry = Path('../../utils/data/')

    file = direcotry / 'bcc_bulk_54_expanded_2_high.xyz'

    for atoms in iread(file.as_posix(), index=slice(1)):
        print(atoms)

        d = atoms2dict(atoms)
        print(d)
        new_atoms = AtomsRow(d).toatoms()
        print(new_atoms)
        print(new_atoms == atoms)

    # d = atoms2dict(atoms)
    # AtomsRow(d).toatoms()
