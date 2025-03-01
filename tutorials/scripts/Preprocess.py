import json
from pathlib import Path

from ase.io import read, write

# import numpy.linalg as la
import matplotlib.pyplot as plt
import numpy as np


class Calculation:
    def __init__(self, *args, **kwargs):
        self.filepath = kwargs.pop("filepath", None)
        self.parameters = kwargs

    def get_data(self, index=-1):
        return read(str(self.filepath))

    @classmethod
    def from_path(cls, path: Path):
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
        # E_gb = calc.parameters.get('E_gb', None)
        #
        # if E_gb is None:
        #     print(calc.filepath)
        #     print(calc.parameters['converged'])
        # else:

        # energy = 16.02 / (2 * calc.parameters['A'] ) * \
        #     (E_gb - potential_energy_per_atom * calc.parameters['n_at'])

        if calc.parameters.get("converged", None):
            # energy = 16.02 / (2 * calc.parameters['A'] ) * \
            #     (calc.parameters.get('E_gb') - potential_energy_per_atom * calc.parameters['n_at'])
            #
            atoms = calc.get_data()
            cell = atoms.get_cell()
            A = cell[0, 0] * cell[1, 1]

            energy = (
                eV
                / Angstrom**2
                / (2 * A)
                * (
                    atoms.get_potential_energy()
                    - potential_energy_per_atom * len(atoms)
                )
            )

            write(calc.filepath.name, atoms)

            # print(energy)
            # print(calc.parameters['converged'])
            # print(data.get_potential_energy())  # data.get_total_energy() == data.get_potential_energy()
            # energies.append(calc.parameters['E_gb'] - data.get_total_energy())
            energies.append(energy)
            angles.append(calc.parameters["angle"] * 180.0 / np.pi)
        else:
            print("not converged: ", calc.filepath)

    plt.bar(angles, energies)

    # x_smooth = np.linspace(min(angles), max(angles), 1000, endpoint=True)
    # f = interp1d(angles, energies, kind='cubic')
    # plt.plot(x_smooth, f(x_smooth), '-')

    plt.show()
