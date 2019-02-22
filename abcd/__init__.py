import logging
from abcd import backends

logger = logging.getLogger(__name__)


class ABCD(object):
    """Factory method"""

    def __new__(cls, url, *args, **kwargs):
        if url.startswith('mongodb://'):
            return backends.MongoDatabase(url=url, *args, **kwargs)

        elif url.startswith('http://') or url.startswith('https://'):
            raise NotImplementedError('http not yet supported! soon...')
        else:
            raise NotImplementedError(f'Unable to recognise the type of connection. (url: {url})')


if __name__ == '__main__':
    from ase.io import iread

    logging.basicConfig(level=logging.INFO)

    for atoms in iread('../utils/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
        # Hack to fix the representation of forces
        atoms.calc.results['forces'] = atoms.arrays['force']

        print(atoms)

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
