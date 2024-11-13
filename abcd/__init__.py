from enum import Enum
import logging
from urllib import parse

logger = logging.getLogger(__name__)


class ConnectionType(Enum):
    mongodb = 1
    http = 2


class ABCD:
    @classmethod
    def from_config(cls, config):
        # Factory method
        url = config["url"]
        return ABCD.from_url(url)

    @classmethod
    def from_url(cls, url, **kwargs):
        # Factory method
        r = parse.urlparse(url)
        logger.info(r)

        if r.scheme == "mongodb":
            conn_settings = {
                "host": r.hostname,
                "port": r.port,
                "username": r.username,
                "password": r.password,
                "authSource": "admin",
            }

            db = r.path.split("/")[1] if r.path else None
            db = db if db else "abcd"

            from abcd.backends.atoms_pymongo import MongoDatabase

            return MongoDatabase(db_name=db, **conn_settings, **kwargs)
        if r.scheme == "mongodb+srv":
            db = r.path.split("/")[1] if r.path else None
            db = db if db else "abcd"
            from abcd.backends.atoms_pymongo import MongoDatabase

            return MongoDatabase(db_name=db, host=r.geturl(), uri_mode=True, **kwargs)
        if r.scheme == "http" or r.scheme == "https":
            raise NotImplementedError("http not yet supported! soon...")
        if r.scheme == "ssh":
            raise NotImplementedError("ssh not yet supported! soon...")
        raise NotImplementedError(
            f"Unable to recognise the type of connection. (url: {url})"
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # url = 'mongodb://mongoadmin:secret@localhost:27017'
    url = "mongodb://mongoadmin:secret@localhost:27017/abcd_new"
    abcd = ABCD.from_url(url)
    abcd.print_info()

    # from ase.io import iread
    # for atoms in iread('../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(1)):
    #     # Hack to fix the representation of forces
    #     atoms.calc.results['forces'] = atoms.arrays['force']
    #
    #     abcd.push(atoms)
    #     print(atoms)
