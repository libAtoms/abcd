import hashlib
import logging
import json
import requests

from typing import Callable, Iterable, Union, Optional, List

# from ase import Atoms
import ase
from ase.calculators.singlepoint import SinglePointCalculator

logger = logging.getLogger(__name__)


class Atoms(ase.Atoms):

    @classmethod
    def from_dict(cls, data):
        return cls(numbers=data['numbers'], positions=data['positions'])


class ABCD(object):
    """client/local interface"""

    def __init__(self, url='http://localhost'):
        #  todo: authentication(user, token), https
        self.url = url

    def push(self, atoms: ase.Atoms):
        # todo: list of Atoms, metadata(user, project, tags)
        message = json.dumps(atoms)
        message_hash = hashlib.md5(message.encode('utf-8')).hexdigest()
        logger.info(message)
        return message_hash

    def pull(self, query=None, properties=None):
        # atoms = json.loads(message)
        atoms = None
        return atoms

    def query(self, query_string):
        pass

    def search(self, query_string: str) -> List[str]:
        results = requests.get(self.url + '/calculations/').json()
        return results

    def get_atoms(self, id: str) -> Atoms:
        data = requests.get(self.url + f'/calculation/{id}').json()
        atoms = Atoms.from_dict(data)
        return atoms

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


if __name__ == '__main__':
    from ase.io import iread

    logging.basicConfig(level=logging.INFO)

    for atoms in iread('../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
        # Hack to fix the representation of forces
        atoms.calc.results['forces'] = atoms.arrays['force']

    print(atoms)

    # todo: query_str2_query_dict

    db = ABCD(url='http://localhost:5000/api')

    query = {
        'elements': ['Cu']
    }
    results = db.search(query)

    # Fetch all results (returns with an Atoms object
    for id in results:
        atoms = db.get_atoms(id)
        print(atoms)

        local_db = [db.get_atoms(id) for id in results]
    #
    # context manager for the database
    with ABCD(url='http://localhost:5000/api') as db:
        results = db.search('formula=Fe3O1;elements=[Fe,*];n_atoms=10,pbc;metadata.collection=iron')
        local_db = [db.get_atoms(id) for id in results]


    # connection = ABCD()
    # for atoms in traj:
    #     hash_value = connection.push(atoms)
    #     print(hash_value)
