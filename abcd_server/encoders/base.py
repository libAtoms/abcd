from abc import ABCMeta, abstractmethod


class BaseEncoder(object, metaclass=ABCMeta):
    """Abstract class for the visitor pattern"""
    default_properties = []

    # @abstractmethod
    def __init__(self):
        pass

    def __enter__(self):
        """support with statement and error handling in python"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def encode(self, atoms):
        """main function"""
        return self.visit_atoms(atoms)

    @abstractmethod
    def visit_atoms(self, atoms):
        pass

    def visit_numbers(self, atoms):
        return atoms.numbers

    def visit_cell(self, atoms):
        return atoms.cell

    def visit_pbc(self, atoms):
        return atoms.pbc

    def visit_positions(self, atoms):
        return atoms.positions

    def visit_forces(self, atoms):
        return atoms.get_forces()

    def visit_energy(self, atoms):
        return atoms.get_potential_energy()

