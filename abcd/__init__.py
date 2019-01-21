import hashlib
import logging
import json
import base64
import numpy as np

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

# import ase

logger = logging.getLogger(__name__)


class ABCD(object):
    """client/local interface"""

    def __init__(self, url='http://localhost'):
        #  todo: authentication(user, token), https
        self.url = url

    def push(self, atoms):
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


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    from ase.io import read

    # from ase.io import write
    # from ase.io.jsonio import read_json
    # write(atoms, '-', format='json')

    traj = read('../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None))

    for atoms in traj:
        # Hack to fix the representation of forces
        atoms.calc.results['forces'] = atoms.arrays['force']

    print(traj)

    atoms = traj[0]
    print(atoms)

    #
    # connection = ABCD()
    # for atoms in traj:
    #     hash_value = connection.push(atoms)
    #     print(hash_value)
