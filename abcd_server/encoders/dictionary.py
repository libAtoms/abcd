from .base import BaseEncoder
import numpy as np


class DictEncoder(BaseEncoder):
    default_properties = ['numbers', 'positions']

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
