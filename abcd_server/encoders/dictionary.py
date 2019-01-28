import ase
import numpy as np

from abcd_server.encoders.base import BaseEncoder, BaseEncoderNew


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


class DictEncoder(BaseEncoderNew):
    # default_properties = ['numbers', 'positions']

    def __init__(self):
        super().__init__()

    def encode(self, atoms: ase.Atoms) -> dict:
        """ASE's original implementation"""
        arrays = atoms.arrays.copy()

        dct = {
            'cell': atoms.cell.tolist(),
            'pbc': atoms.pbc.tolist(),
            'numbers': arrays.pop('numbers').tolist(),
            'positions': arrays.pop('positions').tolist(),
            'arrays': {},
            'info': {},
            'results': {},
            'constraints': [],
        }

        for key, value in arrays.items():

            if isinstance(value, np.ndarray):
                dct['arrays'][key] = value.tolist()
                continue

            dct[key] = value

        for key, value in atoms.info.items():

            if isinstance(value, np.ndarray):
                dct['info'][key] = value.tolist()
                continue

            dct['info'][key] = value

        if atoms.calc is not None:
            dct['results']['calculator_name'] = atoms.calc.name.lower(),
            dct['results']['calculator_parameters'] = atoms.calc.todict()

            for key, value in atoms.calc.results.items():

                if isinstance(value, np.ndarray):
                    dct['results'][key] = value.tolist()
                    continue

                dct['results'][key] = value

        # if atoms.constraints:
        #     dct['constraints'] = [c.todict() for c in atoms.constraints]

        return dct

    def encode_many(self, traj):
        for atoms in traj:
            yield self.encode(atoms)

#
# class DictEncoder(BaseEncoder):
#     # default_properties = ['numbers', 'positions']
#
#     def __init__(self):
#         super().__init__()
#
#     def visit_atoms(self, atoms):
#         """walk through on the structure of an Atoms object """
#         data = {
#             'numbers': self.visit_numbers(atoms),
#             'cell': self.visit_cell(atoms),
#             'pbc': self.visit_pbc(atoms),
#             'positions': self.visit_positions(atoms),
#             'forces': self.visit_forces(atoms),
#             'energy': self.visit_energy(atoms),
#             'metadata': self.visit_metadata(atoms)
#         }
#         return data
#
#     def visit_metadata(self, atoms):
#         # metadata = {}
#         # for key, val in atoms.arrays.items():
#         #     if key not in self.default_properties:
#         #         metadata[key] = self.convert(val)
#         return {key: self.convert(val) for key, val in atoms.arrays.items() if key not in self.default_properties}
#
#     def visit_numbers(self, atoms):
#         return super().visit_numbers(atoms).tolist()
#
#     def visit_cell(self, atoms):
#         return super().visit_cell(atoms).tolist()
#
#     def visit_pbc(self, atoms):
#         return super().visit_pbc(atoms).tolist()
#
#     def visit_positions(self, atoms):
#         return super().visit_positions(atoms).tolist()
#
#     def visit_forces(self, atoms):
#         return super().visit_forces(atoms).tolist()
#
#     def convert(self, data):
#         if isinstance(data, np.ndarray):
#             return data.tolist()
#         if isinstance(data, dict):
#             return {key: self.convert(val) for key, val in data.items()}
#         return data


# TODO: unique id should be a hash value


if __name__ == '__main__':
    from pathlib import Path
    from ase.io import iread

    from ase.db.row import AtomsRow, atoms2dict

    import bson
    from bson.codec_options import CodecOptions

    direcotry = Path('../../utils/data/')

    file = direcotry / 'bcc_bulk_54_expanded_2_high.xyz'

    for atoms in iread(file.as_posix(), index=slice(1)):
        print(atoms)

        with DictEncoder() as encoder:
            d = encoder.encode(atoms)

        print(d)

        new_atoms = AtomsRow(d).toatoms()
        print(new_atoms)
        print(new_atoms == atoms)

    d = encoder.encode(atoms)
    print(bson.BSON.encode(d))
