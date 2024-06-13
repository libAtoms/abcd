from pathlib import Path
from pprint import pprint
import json
from ase.io import read, write
from ase.geometry import crystal_structure_from_cell
import numpy as np

# import numpy.linalg as la

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class Calculation(object):
    def __init__(self, *args, **kwargs):
        self.filepath = kwargs.pop("filepath", None)
        self.parameters = kwargs

    def get_data(self, index=-1):
        return read(str(self.filepath), index=index)

    @classmethod
    def from_path(cls, path: Path, index=-1):
        with (path / "gb.json").open() as data_file:
            gb_data = json.load(data_file)

        with (path / "subgb.json").open() as data_file:
            subgb_data = json.load(data_file)

        # print(gb_data['angle'])

        filename = subgb_data["name"] + "_traj.xyz"
        filepath = (path / filename).resolve()
        # configuration = read(str((path / filename).resolve()), index=index)
        # # gb = read(str((path / filename).resolve()), index=-1)

        # print('{:=^60}'.format(' '+str(path)+' '))
        #
        # print('{:-^40}'.format(' gb.json '))
        # pprint(gb_data)
        #
        # print('{:-^40}'.format(' subgb.json '))
        # pprint(subgb_data)

        # print('{:-^40}'.format(' initial '))
        #
        # force_norm = np.linalg.norm(gb_initial.arrays['force'], axis=1)
        # force_mean = np.mean(force_norm)
        # force_std = np.std(force_norm)
        #
        # print('Force mean: {:f}, std: {:f}'.format(force_mean, force_std))
        # pprint(gb_initial.calc.results)
        #
        # print('{:-^40}'.format(' final '))
        #
        # force_norm = np.linalg.norm(gb_final.arrays['force'], axis=1)
        # force_mean = np.mean(force_norm)
        # force_std = np.std(force_norm)
        #
        # print('Force mean: {:f}, std: {:f}'.format(force_mean, force_std))
        # pprint(gb_final.calc.results)

        return cls(**{**gb_data, **subgb_data}, filepath=filepath)


if __name__ == "__main__":

    # Read grain boundary database
    dirpath = Path("../GB_alphaFe_001")

    calculations = {
        "tilt": [
            Calculation.from_path(calc_dir)
            for calc_dir in (dirpath / "tilt").iterdir()
            if calc_dir.is_dir()
        ],
        "twist": [
            Calculation.from_path(calc_dir)
            for calc_dir in (dirpath / "twist").iterdir()
            if calc_dir.is_dir()
        ],
    }

    # potential energy of the perfect crystal according to a specific potential
    potential_energy_per_atom = -4.01298214176  # alpha-Fe PotBH
    eV = 1.6021766208e-19
    Angstrom = 1.0e-10

    angles, energies = [], []
    for calc in sorted(calculations["tilt"], key=lambda item: item.parameters["angle"]):
        angles.append(calc.parameters["angle"] * 180.0 / np.pi)

        energy = (
            16.02
            / (2 * calc.parameters["A"])
            * (
                calc.parameters["E_gb"]
                - potential_energy_per_atom * calc.parameters["n_at"]
            )
        )

        atoms = calc.get_data()
        cell = atoms.get_cell()
        A = cell[0, 0] * cell[1, 1]

        energy = (
            eV
            / Angstrom**2
            / (2 * A)
            * (atoms.get_total_energy() - potential_energy_per_atom * len(atoms))
        )

        print(energy)
        # print(data.get_potential_energy())  # data.get_total_energy() == data.get_potential_energy()
        # energies.append(calc.parameters['E_gb'] - data.get_total_energy())
        energies.append(energy)

    plt.bar(angles, energies)

    # x_smooth = np.linspace(min(angles), max(angles), 1000, endpoint=True)
    # f = interp1d(angles, energies, kind='cubic')
    # plt.plot(x_smooth, f(x_smooth), '-')

    plt.show()


# 00171501160
# {
#   "A": 128.39243418897394,
#   "gbid": "00171501160",
#   "sigma_csl": 0,
#   "H": 181.47714791126677,
#   "n_at": 514,
#   "orientation_axis": [
#     0,
#     0,
#     1
#   ],
#   "boundary_plane": [
#     1.0,
#     16.0,
#     0.0
#   ],
#   "coincident_sites": 2,
#   "angle": 0.12479104151759456,
#   "type": "symmetric tilt boundary",
#   "zplanes": [
#     45.2801,
#     136.0193
#   ]
# }
# {
#   "A": 1540.7092103178,
#   "E_gb": -98740.56229920377,
#   "gbid": "00171501160_v6bxv2_tv0.4bxv0.1_d1.4z",
#   "name": "00171501160_v6bxv2_tv0.4bxv0.1_d1.4z",
#   "area": 1540.7092103178,
#   "H": 181.47714791,
#   "n_at": 24636,
#   "rbt": [
#     0.4,
#     0.1
#   ],
#   "rcut": 1.4,
#   "converged": true,
#   "param_file": "PotBH.xml",
#   "E_gb_init": -97920.0696568818
# }

# PotBH
#
# Min. Energy: 0.6408 J/m^{2}
# Min. Energy Structure
# Microscopic Degrees of Freedom (x,y, Cutoff Criterion)
# (0.4, 0.1, 1.4)
# (0.4, 0.1, 1.5)
# Max. Energy:1.5965 J/m^{2}
# Max. Energy Structure
# Microscopic Degrees of Freedom (x,y, Cutoff Criterion)
# (0.0, 0.1, 2.1)

#
#
# atoms_mongoengine.py
# gb_ener = 16.02*((sub_dict['E_gb']-(-4.2731*sub_dict['n_at']))/(2*sub_dict['A']))
#
# gb_ener = 16.02*((gb_dict['E_gb']-(ener_per_atom[param_file]*float(gb_dict['n_at'])))/(2*gb_dict['A']))
#
# calc_gap.py
# subgbs = [(16.02*(subgb['E_gb']-float(subgb['n_at']*ener_per_atom['PotBH.xml']))/(2.0*subgb['area']), subgb) for subgb in subgbs]
#
# calc_energy.py
# def calc_e_gb(at, E_bulk):
#   cell = at.get_cell()
#   A    = cell[0,0]*cell[1,1]
#   E_gb = (at.get_potential_energy()-(at.n*(E_bulk)))/(2.*A)
#   print A
#   print at.get_potential_energy()
#   print E_gb, 'eV/A^2'
#   E_gb = 16.02*(at.get_potential_energy()-(at.n*(E_bulk)))/(2.*A)
#   print E_gb, 'J/m^2'
#   return E_gb
#
# relax.py
# E_gb_init   = grain.get_potential_energy()
# E_gb    = grain.get_potential_energy()
#
#
# ener_per_atom = -4.01298214176
