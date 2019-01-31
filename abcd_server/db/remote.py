import logging
import requests
from os import linesep
from typing import List

import ase
from abcd_server.db.base import Database
from abcd_server.formats import DictEncoder

logger = logging.getLogger(__name__)


class Atoms(ase.Atoms):

    @classmethod
    def from_dict(cls, data):
        return cls(numbers=data['numbers'], positions=data['positions'])


class RemoteDatabase(Database):
    """client/local interface"""

    def __init__(self, url='http://localhost'):
        super().__init__()

        self.url = url

    def push(self, atoms: ase.Atoms):
        # todo: list of Atoms, metadata(user, project, tags)

        message = requests.put(self.url + '/calculation', json=DictEncoder().encode(atoms))
        # message = json.dumps(atoms)
        # message_hash = hashlib.md5(message.encode('utf-8')).hexdigest()
        logger.info(message)
        return message.json()

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        raise NotImplementedError

    def query(self, query_string):
        pass

    def search(self, query_string: str) -> List[str]:
        results = requests.get(self.url + '/calculation').json()
        return results

    def get_atoms(self, id: str) -> Atoms:
        data = requests.get(self.url + f'/calculation/{id}').json()
        atoms = Atoms.from_dict(data)
        return atoms

    def __repr__(self):
        return f'ABCD(type={self.connection_type.name}, url={self.url}, ...)'

    def _repr_html_(self):
        """jupyter notebook representation"""
        return '<b>ABCD database</b>'

    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join([
            '{:=^50}'.format(' ABCD Database '),
            '{:>10}: {}'.format('type', 'remote (http/https)'),
            linesep.join('{:>10}: {}'.format(k, v) for k, v in self.db.info().items())
        ])

        print(out)

    def destroy(self):
        self.db.destroy()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
