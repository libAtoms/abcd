import logging
from abcd import backends

from urllib import parse

logger = logging.getLogger(__name__)


class ABCD(object):
    """Factory method"""

    def __new__(cls, url, **kwargs):
        # Factory method

        r = parse.urlparse(url)
        logger.info(r)

        if r.scheme == 'mongodb':

            conn_settings = {
                'host': r.hostname,
                'port': r.port,
                'username': r.username,
                'password': r.password,
                'authentication_source': 'admin',
            }

            db = r.path.split('/')[1] if r.path else None
            db = db if db else 'abcd'

            return backends.MongoDatabase(db=db, **conn_settings, **kwargs)

        elif r.scheme == 'http' or r.scheme == 'https':
            raise NotImplementedError('http not yet supported! soon...')
        else:
            raise NotImplementedError('Unable to recognise the type of connection. (url: {})'.format(url))

    @classmethod
    def from_config(cls, config):
        url = config['url']
        return cls(url)


if __name__ == '__main__':
    from ase.io import iread

    logging.basicConfig(level=logging.INFO)

    url = 'mongodb://mongoadmin:secret@localhost:27017'
    abcd = ABCD(url, collection='atoms')
    abcd.print_info()

    for atoms in iread('../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
        # Hack to fix the representation of forces
        atoms.calc.results['forces'] = atoms.arrays['force']

        abcd.push(atoms)
        print(atoms)
