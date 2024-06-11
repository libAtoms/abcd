from __future__ import annotations

from collections.abc import Iterator
from datetime import datetime
from typing import Iterable, Optional, Union
import logging
from os import linesep
from pathlib import Path

from ase import Atoms
from ase.io import iread
from opensearchpy import (
    OpenSearch,
    helpers,
    AuthenticationException,
    ConnectionTimeout,
    RequestError,
)

from abcd.backends import utils
from abcd.database import AbstractABCD
import abcd.errors
from abcd.model import AbstractModel
from abcd.parsers import extras
from abcd.queryset import AbstractQuerySet


logger = logging.getLogger(__name__)

map_types = {
    bool: "bool",
    float: "float",
    int: "int",
    str: "str",
    datetime: "date",
    dict: "dict",
}


class OpenSearchQuery(AbstractQuerySet):
    """Class to parse and build queries for OpenSearch."""

    def __call__(self, query: Optional[Union[dict, str, list]]) -> Optional[dict]:
        """
        Parses and builds queries for OpenSearch.

        Parameters
        ----------
        query: Optional[Union[dict, str, list]]
            Query to be parsed for OpenSearch. If passed as a dictionary, the query is
            left unchanged. If passed a string or list, the query is treated as a query
            string, based on Lucene query syntax.

        Returns
        -------
        Optional[dict]
            The parsed query for OpenSearch.
        """
        if not query:
            query = self.get_default_query()

        if isinstance(query, str):
            return self.build_query_string(query)
        if isinstance(query, list):
            if len(query) == 0:
                return None
            if query[0] is None:
                return None
            separator = " AND "
            joined_query = separator.join(query)
            return self.build_query_string(joined_query)

        logger.info("parsed query: %s", query)
        return query if query else None

    @staticmethod
    def build_query_string(query: str) -> dict:
        """
        Build query_string (Lucene syntax) query.

        Parameters
        ----------
        query : str
            Query with Lucene syntax.

        Returns
        -------
        dict
            Parsed query for query_string query.
        """
        return {"query_string": {"query": query}}

    @staticmethod
    def get_default_query() -> dict:
        """
        Defines a default OpenSearch query. Currently, matches all documents.

        Returns
        -------
        The default query for OpenSearch.
        """
        return {"match_all": {}}


class AtomsModel(AbstractModel):
    """
    Class to interface between Atoms data and OpenSearch.

    Attributes
    ----------
    _client: Optional[OpenSearch]
        OpenSearch client.
    _index_name: Optional[str]
        OpenSearch index name.
    """

    def __init__(
        self,
        client: Optional[OpenSearch] = None,
        index_name: Optional[str] = None,
        dict: Optional[dict] = None,
    ):
        """
        Initialises class.

        Parameters
        ----------
        client: Optional[OpenSearch]
            OpenSearch client. Default is `None`.
        index_name: Optional[str]
            OpenSearch index name. Default is `None`.
        dict: Optional[dict]
            Dictionary of atoms data. Default is `None`.
        """
        super().__init__(dict)

        self._client = client
        self._index_name = index_name

    @classmethod
    def from_atoms(
        cls,
        client: OpenSearch,
        index_name: str,
        atoms: Atoms,
        extra_info: Optional[dict] = None,
        store_calc: bool = True,
    ) -> AtomsModel:
        """
        Reads and prepares atoms data and extra information for OpenSearch.

        Parameters
        ----------
        client: OpenSearch
            OpenSearch client.
        index_name: str
            OpenSearch index name.
        atoms: Atoms
            Atoms data to be stored.
        extra_info: Optional[dict]
            Extra information to store in the document with the atoms data.
            Default is `None`.
        store_calc: bool, optional
            Whether to store data from the calculator attached to atoms.
            Default is `True`.

        Returns
        -------
        Data from atoms and extra information to be saved in OpenSearch.
        """
        obj = super().from_atoms(atoms, extra_info, store_calc)
        obj._client = client
        obj._index_name = index_name
        return obj

    @property
    def _id(self) -> Optional[str]:
        """
        Get the OpenSearch document ID stored in data.

        Returns
        -------
        Optional[str]
            Current document ID.
        """
        return self.get("_id", None)

    def save(self):
        """
        Saves data in OpenSearch. If the data being saved includes a document
        ID, updates the matching document in OpenSearch with the current data.
        """
        body = {}
        body.update(self.data)
        body["derived"] = self.derived
        if self._client is not None:
            if not self._id:
                self._client.index(index=self._index_name, body=body)
            else:
                del body["_id"]
                body = {"doc": body}
                self._client.update(index=self._index_name, id=self._id, body=body)

    def remove(self):
        """
        If current data includes a document ID, deletes the matching document
        OpenSearch.
        """
        if self._client is not None and self._id:
            self._client.delete(index=self._index_name, id=self._id)
            self.clear()


class OpenSearchDatabase(AbstractABCD):
    """
    Wrapper to make OpenSearch operations easy.

    Attributes
    ----------
    client: OpenSearch
        OpenSearch client.
    db_name: str
        Database name.
    index_name: str
        OpenSearch index name.
    parser: OpenSearchQuery
        Query parser and builder for OpenSearch queries.
    """

    def __init__(
        self,
        host: str = "localhost",
        port: int = 9200,
        db_name: str = "abcd",
        index_name: str = "atoms",
        username: str = "admin",
        password: str = "admin",
        **kwargs,
    ):
        """
        Initialises class.

        Parameters
        ----------
        host: str, optional
            Name of OpenSearch host. Default is `localhost`.
        port: int, optional
            OpenSearch port. Default is `9200`.
        db_name: str, optional
            Label for OpenSearch database. Used only when printing information.
            Default is `abcd`.
        index_name: str, optional
            Name of OpenSearch index. Default is `atoms`.
        username: str, optional
            OpenSearch username. Default is `admin`.
        password: str, optional
            OpenSearch password. Default is `admin`.
        """
        super().__init__()

        logger.info((host, port, index_name, username, password, kwargs))

        client_settings = {
            "verify_certs": False,
            "ca_certs": None,
            "use_ssl": True,
            "ssl_assert_hostname": False,
            "ssl_show_warn": False,
        }

        for key in client_settings:
            if key in kwargs:
                client_settings[key] = kwargs[key]

        self.client = OpenSearch(
            hosts=[{"host": host, "port": port}],
            http_auth=(username, password),
            **client_settings,
        )

        try:
            info = self.client.info()
            logger.info("DB info: %s", info)

        except AuthenticationException as err:
            raise abcd.errors.AuthenticationError() from err

        except ConnectionTimeout as err:
            raise abcd.errors.TimeoutError() from err

        self.db = db_name
        self.index_name = index_name
        self.create()
        self.parser = OpenSearchQuery()

    def info(self):
        """
        Gets information from OpenSearch client about the database.

        Returns
        -------
        Dictionary of database information.
        """
        if self.client.transport.hosts is not None:
            host = self.client.transport.hosts[0]["host"]
            port = self.client.transport.hosts[0]["port"]
        else:
            host, port = None, None

        self.refresh()
        return {
            "host": host,
            "port": port,
            "db": self.db,
            "index": self.index_name,
            "number of confs": self.client.count(index=self.index_name)["count"],
            "type": "opensearch",
        }

    def delete(self, query: Optional[Union[dict, str]] = None):
        """
        Deletes documents from the database.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to be deleted. Default is `None`.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)
        body = {"query": query}

        self.client.delete_by_query(
            index=self.index_name,
            body=body,
        )

    def destroy(self):
        """
        Deletes the current index in OpenSearch.
        Ignores errors if the index does not exist.
        """
        self.client.indices.delete(index=self.index_name, ignore=404)

    def create(self):
        """
        Creates a new index in OpenSearch.
        Ignores errors if the index already exists.
        """
        self.client.indices.create(index=self.index_name, ignore=400)

    def refresh(self):
        """
        Refresh index to ensure recent operations performed are available for search.
        """
        self.client.indices.refresh(index=self.index_name)

    def save_bulk(self, actions: Iterable[dict], **kwargs):
        """
        Save a collection of documents in bulk.

        Parameters
        ----------
        actions: Iterable
            Documents to be saved.
        """
        request_timeout = kwargs.get("request_timeout", 30)
        chunk_size = kwargs.get("chunk_size", 500)
        helpers.bulk(
            client=self.client,
            actions=actions,
            index=self.index_name,
            chunk_size=chunk_size,
            request_timeout=request_timeout,
        )

    def push(
        self,
        atoms: Union[Atoms, Iterable],
        extra_info: Optional[Union[dict, str, list]] = None,
        store_calc: bool = True,
        **kwargs,
    ):
        """
        Save data from atoms object(s) to database.

        Parameters
        ----------
        atoms: Union[Atoms, Iterable]
        extra_info: Optional[Union[dict, str, list]]
            Extra information to store in the document with the atoms data.
            Default is `None`.
        store_calc: bool, optional
            Whether to store data from the calculator attached to atoms.
            Default is `True`.
        """
        if extra_info and isinstance(extra_info, str):
            extra_info = extras.parser.parse(extra_info)  # type: ignore
        if extra_info and isinstance(extra_info, list):
            for i, info in enumerate(extra_info):
                if isinstance(info, str):
                    extra_info[i] = extras.parser.parse(info)

        if isinstance(atoms, Atoms):
            data = AtomsModel.from_atoms(
                self.client,
                self.index_name,
                atoms,
                extra_info=extra_info,  # type: ignore
                store_calc=store_calc,
            )
            data.save()

        elif isinstance(atoms, Iterator) or isinstance(atoms, list):
            actions = []
            for i, item in enumerate(atoms):
                if isinstance(extra_info, list):
                    info = extra_info[i]
                else:
                    info = extra_info
                data = AtomsModel.from_atoms(
                    self.client,
                    self.index_name,
                    item,
                    extra_info=info,  # type: ignore
                    store_calc=store_calc,
                )
                actions.append(data.data)
                actions[-1]["derived"] = data.derived
            self.save_bulk(actions, **kwargs)

    def upload(
        self,
        file: Path,
        extra_infos: Union[Iterable, dict] = (),
        store_calc: bool = True,
    ):
        """
        Upload data from a file to the database.

        Parameters
        ----------
        file: Path
            Path to file to be uploaded
        extra_infos: Union[Iterable, dict]
            Extra information to store in the document with the atoms data.
            Default is `()`.
        store_calc: bool, optional
            Whether to store data from the calculator attached to atoms.
            Default is `True`.
        """

        if isinstance(file, str):
            file = Path(file)

        extra_info = dict(map(extras.parser.parse, extra_infos))

        extra_info["filename"] = str(file)

        data = iread(str(file))
        self.push(data, extra_info, store_calc=store_calc)

    def get_items(self, query: Optional[Union[dict, str]] = None) -> Iterator[dict]:
        """
        Get data as a dictionary from documents in the database.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to get data from. Default is `None`.

        Returns
        -------
        Iterator[dict]
            Iterator for dictionary of data.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)
        query = {
            "query": query,
        }

        for hit in helpers.scan(
            self.client,
            index=self.index_name,
            query=query,
        ):
            yield {"_id": hit["_id"], **hit["_source"]}

    def get_atoms(self, query: Optional[Union[dict, str]] = None) -> Iterator[Atoms]:
        """
        Get data as Atoms object from documents in the database.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to get data from. Default is `None`.

        Returns
        -------
        Iterator[Atoms]
            Generator for AtomsModel object of data.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)
        query = {
            "query": query,
        }

        for hit in helpers.scan(
            self.client,
            index=self.index_name,
            query=query,
        ):
            yield AtomsModel(dict=hit["_source"]).to_ase()

    def count(self, query: Optional[Union[dict, str]] = None, timeout=30.0) -> int:
        """
        Counts number of documents in the database.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to be counted. Default is `None`.
        timeout: float
            Timeout for request in seconds.

        Returns
        -------
        Count of number of documents.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)
        body = {
            "query": query,
        }

        return self.client.count(index=self.index_name, body=body, timeout=timeout)[
            "count"
        ]

    def _get_props_from_source(
        self,
        names: Union[str, list[str]],
        query: Optional[Union[dict, str]] = None,
    ) -> dict:
        """
        Gets all values of specified properties using the original data from _source.

        Parameters
        ----------
        names: Union[str, list[str]]
            Name or list of names of properties to return.
        query: Optional[Union[dict, str]]
            Query to filter documents to get properties from. Default is `None`.

        Returns
        -------
        dict
            Dictionary of lists of values for the specified properties.
        """
        props = {}
        hits = [
            dict(hit["_source"].items())
            for hit in helpers.scan(
                self.client,
                index=self.index_name,
                query=query,
                stored_fields=names,
                _source=names,
            )
            if "_source" in hit and all(name in hit["_source"] for name in names)
        ]
        for name in names:
            props[name] = [hit[name] for hit in hits]
        return props

    def property(
        self,
        names: Union[str, list[str]],
        allow_flatten: bool = True,
        query: Optional[Union[dict, str]] = None,
    ) -> Union[dict, list]:
        """
        Gets all values of specified properties for matching documents in the database.

        Parameters
        ----------
        names: Union[str, list[str]]
            Name or list of names of properties to return.
        allow_flatten: bool = True
            Whether to allow arrays to be returned flattened. There is no guarantee
            for the order of returned values. Default is `True`.
        query: Optional[Union[dict, str]]
            Query to filter documents to get properties from. Default is `None`.

        Returns
        -------
        Union[dict, list]
            Dictionary of lists of values for the specified properties, or list
            if only one property is given.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)
        query = {
            "query": query,
        }

        if isinstance(names, str):
            names = [names]
        names = [format(name) for name in names]

        # Try to use docvalue_fields to avoid loading entire document
        # But not all datatypes supported by default
        if allow_flatten:
            props = {}
            try:
                hits = [
                    dict(hit["fields"].items())
                    for hit in helpers.scan(
                        self.client,
                        index=self.index_name,
                        query=query,
                        _source=False,
                        stored_fields="_none_",
                        docvalue_fields=names,
                    )
                    if "fields" in hit and all(name in hit["fields"] for name in names)
                ]
                for name in names:
                    props[name] = [
                        hit[name][0] if len(hit[name]) == 1 else hit[name]
                        for hit in hits
                    ]

            except RequestError:
                props = self._get_props_from_source(names, query)

        # Use _source to ensure arrays are not flattened
        else:
            props = self._get_props_from_source(names, query)

        if len(names) == 1:
            return props[names[0]]
        return props

    def count_property(self, name, query: Optional[Union[dict, str]] = None) -> dict:
        """
        Counts values of a specified property for matching documents in the
        database. This method much faster than performing a Count on the list
        returned by self.property, so this method should be used preferentially.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to count properties from. Default is `None`.

        Returns
        -------
        Dictionary of values and counts for the specified property for all
        matching documents.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)

        body = {
            "size": 0,
            "query": query,
            "aggs": {
                format(name): {
                    "terms": {
                        "field": format(name),
                        "size": 10000,  # Use composite for all results?
                    },
                },
            },
        }

        prop = {}

        for val in self.client.search(index=self.index_name, body=body)["aggregations"][
            format(name)
        ]["buckets"]:
            prop[val["key"]] = val["doc_count"]

        return prop

    def properties(self, query: Optional[Union[dict, str]] = None) -> dict:
        """
        Gets lists of all properties from matching documents, separated into
        info, derived, and array properties.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to get properties from. Default is `None`.

        Returns
        -------
        Dictionary of properties, with keys corresponding to info, derived,
        and arrays of properties, and values corresponding to a list of
        the properties of that type.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)

        properties = {}

        for prop in self.client.indices.get_mapping(index=self.index_name)[
            self.index_name
        ]["mappings"]["properties"].keys():
            body = {
                "size": 0,
                "query": query,
                "aggs": {
                    "info_keys": {
                        "filter": {"term": {"derived.info_keys.keyword": prop}},
                    },
                    "derived_keys": {
                        "filter": {"term": {"derived.derived_keys.keyword": prop}},
                    },
                    "arrays_keys": {
                        "filter": {"term": {"derived.arrays_keys.keyword": prop}},
                    },
                },
            }

            res = self.client.search(
                index=self.index_name,
                body=body,
            )

            for label in ("info_keys", "derived_keys", "arrays_keys"):
                count = res["aggregations"][label]["doc_count"]
                if count > 0:
                    key = label.split("_", maxsplit=1)[0]
                    if key in properties:
                        properties[key].append(prop)
                    else:
                        properties[key] = [prop]

        return properties

    def get_type_of_property(self, prop: str, category: str) -> str:
        """
        Gets type of a property, given its category.

        Parameters
        ----------
        prop: str
            Name of the property.
        catagory: str
            Name of property's category. Current options are `info`, `derived`,
            and `arrays`.

        Returns
        -------
        Type of the property.
        """
        atoms = self.client.search(
            index=self.index_name,
            body={"size": 1, "query": {"exists": {"field": prop}}},
        )

        data = atoms["hits"]["hits"][0]["_source"][prop]

        if category == "arrays":
            if isinstance(data[0], list):
                return "array({}, N x {})".format(
                    map_types[type(data[0][0])], len(data[0])
                )
            return "vector({}, N)".format(map_types[type(data[0])])

        if isinstance(data, list):
            if isinstance(data[0], list):
                if isinstance(data[0][0], list):
                    return "list(list(...)"
                return "array({})".format(map_types[type(data[0][0])])
            return "vector({})".format(map_types[type(data[0])])
        return "scalar({})".format(map_types[type(data)])

    def count_properties(self, query: Optional[Union[dict, str]] = None) -> dict:
        """
        Counts all properties from matching documents.

        Parameters
        ----------
        query: Optional[Union[dict, str]]
            Query to filter documents to count properties from. Default is `None`.

        Returns
        -------
        Dictionary of properties, with keys property names, and values
        corresponding to their counts, categories and data types.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)
        properties = {}

        try:
            keys = self.client.indices.get_mapping(index=self.index_name)[
                self.index_name
            ]["mappings"]["properties"].keys()
        except KeyError:
            return properties

        for key in keys:
            body = {
                "size": 0,
                "query": query,
                "aggs": {
                    "info_keys": {
                        "filter": {"term": {"derived.info_keys.keyword": key}},
                    },
                    "derived_keys": {
                        "filter": {"term": {"derived.derived_keys.keyword": key}},
                    },
                    "arrays_keys": {
                        "filter": {"term": {"derived.arrays_keys.keyword": key}},
                    },
                },
            }

            res = self.client.search(
                index=self.index_name,
                body=body,
            )

            for label in ("info_keys", "derived_keys", "arrays_keys"):
                count = res["aggregations"][label]["doc_count"]
                if count > 0:
                    properties[key] = {
                        "count": count,
                        "category": label.split("_", maxsplit=1)[0],
                        "dtype": self.get_type_of_property(
                            key, label.split("_", maxsplit=1)[0]
                        ),
                    }

        return properties

    def add_property(self, data: dict, query: Optional[Union[dict, str]] = None):
        """
        Adds properties to matching documents.

        Parameters
        ----------
        data: dict
            Property key-value pairs to be added to matching documents.
        query: Optional[Union[dict, str]]
            Query to filter documents to add properties to. Default is `None`.
        """
        query = self.parser(query)
        logger.info("add: data=%s, query=%s", data, query)

        script_txt = "ctx._source.derived.info_keys.addAll(params.keys);"
        for key, val in data.items():
            script_txt += f"ctx._source.{key} = '{val}';"

        body = {
            "script": {
                "source": script_txt,
                "lang": "painless",
                "params": {"keys": list(data.keys())},
            },
            "query": query,
        }

        self.client.update_by_query(
            index=self.index_name,
            body=body,
        )

    def rename_property(
        self, name: str, new_name: str, query: Optional[Union[dict, str]] = None
    ):
        """
        Renames property for all matching documents.

        Parameters
        ----------
        name: str
            Current name of property to be renamed.
        new_name: str
            New name of property to be renamed.
        query: Optional[Union[dict, str]]
            Query to filter documents to rename property. Default is `None`.
        """
        query = self.parser(query)
        logger.info("rename: query=%s, old=%s, new=%s", query, name, new_name)

        script_txt = "if (!ctx._source.containsKey(params.new_name)) { "
        script_txt += (
            f"ctx._source.{new_name} = ctx._source.{name};"
            " ctx._source.remove(params.name);"
            " for (int i=0; i<ctx._source.derived.info_keys.length; i++) {"
            " if (ctx._source.derived.info_keys[i] == params.name) { "
            " ctx._source.derived.info_keys[i] = params.new_name;}}}"
        )

        body = {
            "script": {
                "source": script_txt,
                "lang": "painless",
                "params": {"name": name, "new_name": new_name},
            },
            "query": query,
        }

        self.client.update_by_query(index=self.index_name, body=body)

    def delete_property(self, name: str, query: Optional[Union[dict, str]] = None):
        """
        Deletes property from all matching documents.

        Parameters
        ----------
        name: str
            Name of property to be deleted from documents.
        query: Optional[Union[dict, str]]
            Query to filter documents to have property deleted. Default is `None`.
        """
        query = self.parser(query)
        logger.info("delete: query=%s, porperty=%s", name, query)

        script_txt = f"if (ctx._source.containsKey('{name}')) {{ "
        script_txt += "ctx._source.remove(params.name);"
        script_txt += "for (int i=0; i<ctx._source.derived.info_keys.length; i++) {"
        script_txt += "if (ctx._source.derived.info_keys[i] == params.name) { "
        script_txt += "ctx._source.derived.info_keys.remove(i);}}}"

        body = {
            "script": {
                "source": script_txt,
                "lang": "painless",
                "params": {
                    "name": name,
                },
            },
            "query": query,
        }

        self.client.update_by_query(index=self.index_name, body=body)

    def hist(
        self, name: str, query: Optional[Union[dict, str]] = None, **kwargs
    ) -> Optional[dict]:
        """
        Calculate histogram statistics for a property from all matching documents.

        Parameters
        ----------
        name: str
            Name of property.
        query: Optional[Union[dict, str]]
            Query to filter documents. Default is `None`.

        Returns
        -------
        Optional[dict]
            Dictionary containing histogram statistics, including the number of
            bins, edges, counts, min, max, and standard deviation.
        """
        query = self.parser(query)
        logger.info("parsed query: %s", query)

        data = self.property(name, query)
        return utils.histogram(name, data, **kwargs)

    def __repr__(self):
        """
        OpenSearch class representation.

        Returns
        -------
        String for OpenSearch class representation, containing the connected
        database host, port, and index name.
        """
        if self.client.transport.hosts is not None:
            host = self.client.transport.hosts[0]["host"]
            port = self.client.transport.hosts[0]["port"]
        else:
            host, port = None, None

        return (
            f"{self.__class__.__name__}("
            f"url={host}:{port}, "
            f"index={self.index_name}) "
        )

    def _repr_html_(self):
        """
        Jupyter notebook representation of OpenSearch class.

        Returns
        -------
        String for HTML representation.
        """
        return "<b>ABCD OpenSearch database</b>"

    def print_info(self):
        """
        Show basic information about the connected OpenSearch database.
        """
        out = linesep.join(
            [
                "{:=^50}".format(" ABCD OpenSearch "),
                "{:>10}: {}".format("type", "opensearch"),
                linesep.join("{:>10}: {}".format(k, v) for k, v in self.info().items()),
            ]
        )

        print(out)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


if __name__ == "__main__":
    db = OpenSearchDatabase(username="admin", password="admin")
    print(db.info())
