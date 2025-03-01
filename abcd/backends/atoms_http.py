# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar
# This program is distributed under the MIT License, see LICENSE.md.

import json
import logging
from os import linesep

import ase
import requests

from abcd.backends.abstract import Database

logger = logging.getLogger(__name__)

DictEncoder = json.JSONEncoder


class Atoms(ase.Atoms):
    @classmethod
    def from_dict(cls, data):
        return cls(numbers=data["numbers"], positions=data["positions"])


class HttpDatabase(Database):
    """client/local interface"""

    def __init__(self, url="http://localhost"):
        super().__init__()

        self.url = url

    def push(self, atoms: ase.Atoms):
        # todo: list of Atoms, metadata(user, project, tags)

        message = requests.put(
            self.url + "/calculation", json=DictEncoder().encode(atoms)
        )
        # message = json.dumps(atoms)
        # message_hash = hashlib.md5(message.encode('utf-8')).hexdigest()
        logger.info(message)
        return message.json()

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        raise NotImplementedError

    def query(self, query_string):
        pass

    def search(self, query_string: str) -> list[str]:
        return requests.get(self.url + "/calculation").json()

    def get_atoms(self, id: str) -> Atoms:
        data = requests.get(self.url + f"/calculation/{id}").json()
        return Atoms.from_dict(data)

    def __repr__(self):
        return f"ABCD(type={self.__class__.__name__}, url={self.url}, ...)"

    def _repr_html_(self):
        """jupyter notebook representation"""
        return "<b>ABCD database</b>"

    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join(
            [
                "{:=^50}".format(" ABCD Database "),
                "{:>10}: {}".format("type", "remote (http/https)"),
                linesep.join(f"{k:>10}: {v}" for k, v in self.db.info().items()),
            ]
        )

        print(out)

    def destroy(self):
        self.db.destroy()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


if __name__ == "__main__":
    abcd = HttpDatabase(url="http://localhost:8080/api")
    abcd.print_info()
