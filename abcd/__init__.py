import logging
from urllib import parse
from abc import ABCMeta, abstractmethod
from enum import Enum

logger = logging.getLogger(__name__)


class ConnectionType(Enum):
    mongodb = 1
    http = 2


def from_config(config):
    # Factory method
    url = config['url']
    return from_url(url)


def from_url(url, **kwargs):
    # Factory method
    r = parse.urlparse(url)
    logger.info(r)

    if r.scheme == 'mongodb':

        conn_settings = {
            'host': r.hostname,
            'port': r.port,
            'username': r.username,
            'password': r.password,
            'authSource': 'admin',
        }

        db = r.path.split('/')[1] if r.path else None
        db = db if db else 'abcd'

        from abcd.backends.atoms_pymongo import MongoDatabase
        return MongoDatabase(db_name=db, **conn_settings, **kwargs)

    elif r.scheme == 'http' or r.scheme == 'https':
        raise NotImplementedError('http not yet supported! soon...')
    elif r.scheme == 'ssh':
        raise NotImplementedError('ssh not yet supported! soon...')
    else:
        raise NotImplementedError('Unable to recognise the type of connection. (url: {})'.format(url))


class ABCD(metaclass=ABCMeta):
    """Factory method"""

    @abstractmethod
    def __init__(self):
        pass

    def info(self):
        pass

    def push(self, atoms):
        pass

    def pull(self, query=None, properties=None):
        pass

    def query(self, query_string):
        pass

    def destroy(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __repr__(self):
        pass

    def _repr_html_(self):
        pass

    def print_info(self):
        pass


class QuerySet(metaclass=ABCMeta):
    pass


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # url = 'mongodb://mongoadmin:secret@localhost:27017'
    url = 'mongodb://mongoadmin:secret@localhost:27017/abcd_new'
    abcd = from_url(url)
    abcd.print_info()

    # from ase.io import iread
    # for atoms in iread('../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
    #     # Hack to fix the representation of forces
    #     atoms.calc.results['forces'] = atoms.arrays['force']
    #
    #     abcd.push(atoms)
    #     print(atoms)
