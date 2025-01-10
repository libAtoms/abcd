import datetime
import getpass
import logging
from hashlib import md5
from collections import Counter, UserDict
from ase.calculators.singlepoint import SinglePointCalculator

import numpy as np
from ase import Atoms

logger = logging.getLogger(__name__)


class Hasher(object):
    def __init__(self, method=md5):
        self.method = method()

    def update(self, value):

        if isinstance(value, int):
            self.update(str(value).encode("ascii"))

        elif isinstance(value, str):
            self.update(value.encode("utf-8"))

        elif isinstance(value, float):
            self.update("{:.8e}".format(value).encode("ascii"))

        elif isinstance(value, (tuple, list)):
            for e in value:
                self.update(e)

        elif isinstance(value, (dict, UserDict)):
            keys = value.keys()
            for k in sorted(keys):
                self.update(k.encode("utf-8"))
                self.update(value[k])

        elif isinstance(value, datetime.datetime):
            self.update(str(value))

        elif isinstance(value, bytes):
            self.method.update(value)

        else:
            raise ValueError(
                "The {} type cannot be hashed! (Value: {})", format(type(value), value)
            )

    def __call__(self):
        """Retrieve the digest of the hash."""
        return self.method.hexdigest()


class AbstractModel(UserDict):
    reserved_keys = {
        "n_atoms",
        "cell",
        "pbc",
        "calculator_name",
        "calculator_parameters",
        "derived",
    }

    def __init__(self, dict=None, **kwargs):
        self.arrays_keys = []
        self.info_keys = []
        self.results_keys = []
        self.derived_keys = []

        super().__init__(dict, **kwargs)

    @property
    def derived(self):
        return {
            "arrays_keys": self.arrays_keys,
            "info_keys": self.info_keys,
            "results_keys": self.results_keys,
            "derived_keys": self.derived_keys,
        }

    def __getitem__(self, key):

        if key == "derived":
            return self.derived

        return super().__getitem__(key)

    def __setitem__(self, key, value):

        if key == "derived":
            # raise KeyError('Please do not use "derived" as key because it is protected!')
            # Silent return to avoid raising error in pymongo package
            return

        value = self.convert(value)
        self.update_key_category(key, value)

        super().__setitem__(key, value)

    def convert(self, value):
        # TODO: https://api.mongodb.com/python/current/api/bson/index.html using type_registry

        if isinstance(value, np.int64):
            return int(value)

        return value

    def update_key_category(self, key, value):

        if key == "_id":
            # raise KeyError('Please do not use "derived" as key because it is protected!')
            return

        for category in ("arrays_keys", "info_keys", "results_keys", "derived_keys"):
            if key in self.derived[category]:
                return

        if key in ("positions", "forces"):
            self.derived["arrays_keys"].append(key)
            return

        if key in ("n_atoms", "cell", "pbc"):
            self.derived["info_keys"].append(key)
            return

        # Guess the category based in the shape of the value
        n_atoms = self["n_atoms"]
        if isinstance(value, (np.ndarray, list)) and len(value) == n_atoms:
            self.derived["arrays_keys"].append(key)
        else:
            self.derived["info_keys"].append(key)

    def __delitem__(self, key):
        for category in ("arrays_keys", "info_keys", "results_keys", "derived_keys"):
            if key in self.derived[category]:
                self.derived[category].remove(key)
                break

        super().__delitem__(key)

    def __iter__(self):
        for item in super().__iter__():
            yield item
        yield "derived"

    @classmethod
    def from_atoms(cls, atoms: Atoms, extra_info=None, store_calc=True):
        """Extract data from Atoms info, arrays and results."""
        if not isinstance(atoms, Atoms):
            raise ValueError("atoms must be an ASE Atoms object.")

        reserved_keys = {
            "n_atoms",
            "cell",
            "pbc",
            "calculator_name",
            "calculator_parameters",
            "derived",
            "formula",
        }

        arrays_keys = set(atoms.arrays.keys())
        info_keys = set(atoms.info.keys())
        if store_calc and atoms.calc:
            results_keys = atoms.calc.results.keys() - (arrays_keys | info_keys)
        else:
            results_keys = set()

        all_keys = (reserved_keys, arrays_keys, info_keys, results_keys)
        if len(set.union(*all_keys)) != sum(map(len, all_keys)):
            print(all_keys)
            raise ValueError("All the keys must be unique!")

        item = cls()

        n_atoms = len(atoms)

        data = {
            "n_atoms": n_atoms,
            "cell": atoms.cell.tolist(),
            "pbc": atoms.pbc.tolist(),
            "formula": atoms.get_chemical_formula(),
        }

        info_keys.update(data.keys())

        for key, value in atoms.arrays.items():
            if isinstance(value, np.ndarray):
                data[key] = value.tolist()
            else:
                data[key] = value

        for key, value in atoms.info.items():
            if isinstance(value, np.ndarray):
                data[key] = value.tolist()
            else:
                data[key] = value

        if store_calc and atoms.calc:
            data["calculator_name"] = atoms.calc.__class__.__name__
            data["calculator_parameters"] = atoms.calc.todict()
            info_keys.update({"calculator_name", "calculator_parameters"})

            for key, value in atoms.calc.results.items():
                if isinstance(value, np.ndarray):
                    data[key] = value.tolist()
                else:
                    data[key] = value

        item.arrays_keys = list(arrays_keys)
        item.info_keys = list(info_keys)
        item.results_keys = list(results_keys)

        item.update(data)

        if extra_info:
            item.info_keys.extend(extra_info.keys())
            item.update(extra_info)

        item.pre_save()
        return item

    def to_ase(self):
        arrays_keys = set(self.arrays_keys)
        info_keys = set(self.info_keys)

        cell = self.pop("cell", None)
        pbc = self.pop("pbc", None)
        numbers = self.pop("numbers", None)
        positions = self.pop("positions", None)
        results_keys = self.derived["results_keys"]

        info_keys -= {"cell", "pbc"}
        arrays_keys -= {"numbers", "positions"}

        atoms = Atoms(cell=cell, pbc=pbc, numbers=numbers, positions=positions)

        if "calculator_name" in self:
            # calculator_name = self['info'].pop('calculator_name')
            # atoms.calc = get_calculator(data['results']['calculator_name'])(**params)

            params = self.pop("calculator_parameters", {})
            info_keys -= {"calculator_parameters"}

            atoms.calc = SinglePointCalculator(atoms, **params)
            atoms.calc.results.update((key, self[key]) for key in results_keys)

        atoms.arrays.update((key, np.array(self[key])) for key in arrays_keys)
        atoms.info.update((key, self[key]) for key in info_keys)

        return atoms

    def pre_save(self):
        self.derived_keys = ["elements", "username", "uploaded", "modified"]

        cell = self["cell"]

        if cell:
            volume = abs(np.linalg.det(cell))  # atoms.get_volume()
            self.derived_keys.append("volume")
            self["volume"] = volume

            virial = self.get("virial")
            if virial:
                # pressure P = -1/3 Tr(stress) = -1/3 Tr(virials/volume)
                self.derived_keys.append("pressure")
                self["pressure"] = -1 / 3 * np.trace(virial / volume)

        # 'elements': Counter(atoms.get_chemical_symbols()),
        self["elements"] = Counter(str(element) for element in self["numbers"])

        self["username"] = getpass.getuser()

        if not self.get("uploaded"):
            self["uploaded"] = datetime.datetime.now(datetime.timezone.utc)

        self["modified"] = datetime.datetime.now(datetime.timezone.utc)

        hasher = Hasher()

        for key in ("numbers", "positions", "cell", "pbc"):
            hasher.update(self[key])

        self.derived_keys.append("hash_structure")
        self["hash_structure"] = hasher()

        for key in self.arrays_keys:
            hasher.update(self[key])
        for key in self.info_keys:
            hasher.update(self[key])

        self.derived_keys.append("hash")
        self["hash"] = hasher()


if __name__ == "__main__":
    import io
    from pprint import pprint
    from ase.io import read

    logging.basicConfig(level=logging.INFO)
    # from ase.io import jsonio

    atoms = read("test.xyz", format="xyz", index=0)
    atoms.set_cell([1, 1, 1])

    print(atoms)
    # print(atoms.arrays)
    # print(atoms.info)

    # pprint(AbstractModel.from_atoms(atoms))

    # pprint(jsonio.encode(atoms.arrays))
    # pprint(jsonio.encode(atoms.info))
    # pprint(jsonio.encode(atoms.cell))

    pprint(AbstractModel.from_atoms(atoms))

    h = Hasher()
    h.update(AbstractModel.from_atoms(atoms))
    print(h())

    model = AbstractModel.from_atoms(atoms)
    print(model.to_ase())

    # xyz = io.StringIO(
    #     """
    #     2
    #     Properties=species:S:1:pos:R:3 s="sadf" _vtk_test="t e s t _ s t r" pbc="F F F"
    #     Si       0.00000000       0.00000000       0.00000000
    #     Si       0.00000000       0.00000000       0.00000000
    #
    #     """)
    #
    # atoms = read(xyz, format='extxyz', index=0)
