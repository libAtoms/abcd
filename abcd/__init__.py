import logging
from urllib import parse
from enum import Enum

logger = logging.getLogger(__name__)


class ConnectionType(Enum):
    mongodb = 1
    http = 2
    opensearch = 3


class ABCD(object):
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

        if ConnectionType[r.scheme] is ConnectionType.mongodb:
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

        if ConnectionType[r.scheme] is ConnectionType.opensearch:
            conn_settings = {
                "host": r.hostname,
                "port": r.port,
                "username": r.username,
                "password": r.password,
            }

            db = r.path.split("/")[1] if r.path else None
            db = db if db else "abcd"

            from abcd.backends.atoms_opensearch import OpenSearchDatabase

            return OpenSearchDatabase(db=db, **conn_settings, **kwargs)

        elif r.scheme == "http" or r.scheme == "https":
            raise NotImplementedError("http not yet supported! soon...")
        elif r.scheme == "ssh":
            raise NotImplementedError("ssh not yet supported! soon...")
        else:
            raise NotImplementedError(
                "Unable to recognise the type of connection. (url: {})".format(url)
            )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # url = 'mongodb://mongoadmin:secret@localhost:27017'
    url = "mongodb://mongoadmin:secret@localhost:27017/abcd_new"
    abcd = ABCD.from_url(url)
    abcd.print_info()
