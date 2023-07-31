import logging

from os import linesep

from ase import Atoms

import abcd.errors
from abcd.model import AbstractModel
from abcd.database import AbstractABCD

from opensearchpy import OpenSearch, AuthenticationException, ConnectionTimeout

logger = logging.getLogger(__name__)


class AtomsModel(AbstractModel):
    def __init__(self, client=None, index_name=None, dict=None):
        super().__init__(dict)

        self._client = client
        self._index_name = index_name

    @classmethod
    def from_atoms(cls, client, index_name, atoms: Atoms, extra_info=None, store_calc=True):
        obj = super().from_atoms(atoms, extra_info, store_calc)
        obj._client = client
        obj._index_name = index_name
        return obj

    @property
    def _id(self):
        return self.get("_id", None)


class OpenSearchDatabase(AbstractABCD):
    """Wrapper to make database operations easy"""

    def __init__(
            self,
            host="localhost",
            port=9200,
            index_name="atoms",
            username="admin",
            password="admin",
            **kwargs):

        super().__init__()

        logger.info((host, port, index_name, username, password, kwargs))

        self.client = OpenSearch(
            hosts=[{"host": host, "port": port}],
            http_auth=(username, password),
            verify_certs=False,
            ca_certs=False,
            use_ssl=True,
            ssl_assert_hostname=False,
            ssl_show_warn=False,
        )

        try:
            info = self.client.info()
            logger.info("DB info: {}".format(info))

        except AuthenticationException:
            raise abcd.errors.AuthenticationError()

        except ConnectionTimeout:
            raise abcd.errors.TimeoutError()

        self.index_name = index_name
        self.create()

    def info(self):
        host = self.client.transport.hosts[0]["host"]
        port = self.client.transport.hosts[0]["port"]

        self.client.indices.refresh(index=self.index_name)
        return {
            "host": host,
            "port": port,
            "index": self.index_name,
            "number of confs": self.client.count(index=self.index_name)["count"],
            "type": "opensearch"
        }

    def create(self):
        self.client.indices.create(index=self.index_name, ignore=400)

    def __repr__(self):
        host = self.client.transport.hosts[0]["host"]
        port = self.client.transport.hosts[0]["port"]

        return "{}(".format(self.__class__.__name__) + \
               "url={}:{}, ".format(host, port) + \
               "index={}) ".format(self.index_name)

    def _repr_html_(self):
        """Jupyter notebook representation"""
        return "<b>ABCD OpenSearch database</b>"

    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join(["{:=^50}".format(" ABCD OpenSearch "),
                            "{:>10}: {}".format("type", "opensearch"),
                            linesep.join("{:>10}: {}".format(k, v) for k, v in self.info().items())])

        print(out)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


if __name__ == "__main__":
    db = OpenSearchDatabase(username="admin", password="admin")
    print(db.info())
