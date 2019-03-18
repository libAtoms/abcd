import logging
from urllib import parse
from abcd import backends

logger = logging.getLogger(__name__)


class ABCD(object):
    """Factory method"""

    def __new__(cls, url, **kwargs):

        r = parse.urlparse(url)
        print(r)
        if r.scheme == 'mongodb':

            conn_settings = {
                'host': r.hostname,
                'port': r.port,
                'username': r.username,
                'password': r.password,
                'authentication_source': 'admin',
            }

            db = r.path.split('/', 1)[1]
            db = db if db else 'abcd'

            return backends.MongoDatabase(db=db, **conn_settings, **kwargs)

        elif r.scheme == 'http' or r.scheme == 'https':
            raise NotImplementedError('http not yet supported! soon...')
        else:
            raise NotImplementedError(f'Unable to recognise the type of connection. (url: {url})')

    @classmethod
    def from_config(cls, config):
        url = config['url']
        return cls(url)


if __name__ == '__main__':
    from ase.io import iread

    logging.basicConfig(level=logging.INFO)

    url = 'mongodb://2ef35d3635e9dc5a922a6a42:ac6ce72e259f5ddcc8dd5178@localhost:27017/abcd'
    abcd = ABCD(url, collection='atoms')
    abcd.print_info()

    for atoms in iread('../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
        # Hack to fix the representation of forces
        atoms.calc.results['forces'] = atoms.arrays['force']

        print(atoms)

    abcd.query('aa>22 bb')
    # todo: query_str2_query_dict

    # db = ABCD(url='http://localhost:5000/api')
    #
    # query = {
    #     'elements': ['Cu']
    # }
    # results = db.search(query)
    #
    # # Fetch all results (returns with an Atoms object
    # for id in results:
    #     atoms = db.get_atoms(id)
    #     print(atoms)
    #
    #     local_db = [db.get_atoms(id) for id in results]
    # #
    # # context manager for the database
    # with ABCD(url='http://localhost:5000/api') as db:
    #     results = db.search('formula=Fe3O1;elements=[Fe,*];n_atoms=10,pbc;metadata.collection=iron')
    #     local_db = [db.get_atoms(id) for id in results]

    # with ABCD(url='http://localhost:5000/api') as db:
    #     results = db.search('formula=Fe3O1;elements=[Fe,*];n_atoms=10,pbc;metadata.collection=iron')
    #     local_db = [db.get_atoms(id) for id in results]

    # connection = ABCD()
    # for atoms in traj:
    #     hash_value = connection.push(atoms)
    #     print(hash_value)

    # with ABCD(url='http://localhost:5000/api') as db:
    #     db.push(atoms)
