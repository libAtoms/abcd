from __future__ import annotations

import logging

from collections.abc import Generator
from typing import Union, Iterable
from os import linesep
from datetime import datetime
from collections import Counter
from operator import itemgetter
from pathlib import Path

import numpy as np

from ase import Atoms
from ase.io import iread

import abcd.errors
from abcd.model import AbstractModel
from abcd.database import AbstractABCD
from abcd.queryset import AbstractQuerySet
from abcd.parsers import extras

from opensearchpy import OpenSearch, helpers, AuthenticationException, ConnectionTimeout

from luqum.parser import parser
from luqum.elasticsearch import SchemaAnalyzer, ElasticsearchQueryBuilder

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
    """
    Class to parse and build queries for OpenSearch.

    Attributes
    ----------
    query_builder: ElasticsearchQueryBuilder
        Query builder to convert a Tree in an OpenSearch query.
    """

    def __init__(
        self,
        client: Union[OpenSearch, None] = None,
        index_name: Union[str, None] = None,
        analyse_schema: bool = False,
    ):
        """ "
        Initialises class.

        Parameters
        ----------
        client: Union[OpenSearch, None]
            OpenSearch client, used for if analyse_schema is `True` to
            characterise the schema. Default is `None`.
        index_name: Union[str, None]
            Name of OpenSearch index to be analysed, used if analyse_schema
            is `True` to characterise the schema. Default is `None`.
        analyse_schema: bool, optional
            Whether to analyse the schema, as defined by the index_name and client.
            Default is `False`.
        """
        if analyse_schema and client is not None and index_name is not None:
            schema = client.indices.get_mapping()[index_name]
            schema_analizer = SchemaAnalyzer(schema)
            self.query_builder = ElasticsearchQueryBuilder(
                **schema_analizer.query_builder_options()
            )
        else:
            self.query_builder = ElasticsearchQueryBuilder()

    def __call__(self, query: Union[dict, str, None]) -> Union[dict, None]:
        """
        Parses and builds queries from strings using ElasticsearchQueryBuilder.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to be parsed for OpenSearch. If given as a dictionary,
            the query is left unchanged. If given as a string, the
            ElasticsearchQueryBuilder is used to build the query.

        Returns
        -------
        Union[dict, None]
            The parsed query for OpenSearch.
        """
        logger.info("parsed query: {}".format(query))

        if not query:
            query = self.get_default_query()

        if isinstance(query, dict):
            return query
        elif isinstance(query, str):
            tree = parser.parse(query)
            return self.query_builder(tree)

        return query if query else None

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
    _client: Union[OpenSearch, None]
        OpenSearch client.
    _index_name: Union[str, None]
        OpenSearch index name.
    """

    def __init__(
        self,
        client: Union[OpenSearch, None] = None,
        index_name: Union[str, None] = None,
        dict: Union[dict, None] = None,
    ):
        """
        Initialises class.

        Parameters
        ----------
        client: Union[OpenSearch, None]
            OpenSearch client.
        index_name: Union[str, None]
            OpenSearch index name.
        dict: dict
            Dictionary of atoms data.
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
        extra_info: Union[dict, None] = None,
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
        extra_info: Union[dict, None], optional
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
    def _id(self):
        """
        Gets the OpenSearch document ID stored in data.

        Returns
        -------
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
                body.pop("_id", None)
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
    db: str
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
        db: str = "abcd",
        index_name: str = "atoms",
        username: str = "admin",
        password: str = "admin",
        analyse_schema: bool = True,
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
        db: str, optional
            Label for OpenSearch database. Used only when printing information.
            Default is `abcd`.
        index_name: str, optional
            Name of OpenSearch index. Default is `atoms`.
        username: str, optional
            OpenSearch username. Default is `admin`.
        password: str, optional
            OpenSearch password. Default is `admin`.
        analyse_schema: bool, optional
            Whether to analyse the OpenSearch schema when building queries.
            Default is `True`.
        """
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

        self.db = db
        self.index_name = index_name
        self.create()
        self.parser = OpenSearchQuery(self.client, self.index_name, analyse_schema)

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

        self.client.indices.refresh(index=self.index_name)
        return {
            "host": host,
            "port": port,
            "db": self.db,
            "index": self.index_name,
            "number of confs": self.client.count(index=self.index_name)["count"],
            "type": "opensearch",
        }

    def delete(self, query: Union[dict, str, None] = None):
        """
        Deletes documents from the database.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to be deleted. Default is `None`.
        """
        query = self.parser(query)
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

    def save_bulk(self, actions: Iterable):
        """
        Save a collection of documents in bulk.

        Parameters
        ----------
        actions: Iterable
            Documents to be saved.
        """
        helpers.bulk(client=self.client, actions=actions, index=self.index_name)

    def push(
        self,
        atoms: Union[Atoms, Iterable],
        extra_info: Union[dict, str, None] = None,
        store_calc: bool = True,
    ):
        """
        Save data from atoms object(s) to database.

        Parameters
        ----------
        atoms: Union[Atoms, Iterable]
        extra_info: Union[dict, str, None], optional
            Extra information to store in the document with the atoms data.
            Default is `None`.
        store_calc: bool, optional
            Whether to store data from the calculator attached to atoms.
            Default is `True`.
        """
        if extra_info and isinstance(extra_info, str):
            extra_info = extras.parser.parse(extra_info)  # type: ignore

        if isinstance(atoms, Atoms):
            data = AtomsModel.from_atoms(
                self.client,
                self.index_name,
                atoms,
                extra_info=extra_info,  # type: ignore
                store_calc=store_calc,
            )
            data.save()

        elif isinstance(atoms, Generator) or isinstance(atoms, list):
            actions = []
            for item in atoms:
                data = AtomsModel.from_atoms(
                    self.client,
                    self.index_name,
                    item,
                    extra_info=extra_info,  # type: ignore
                    store_calc=store_calc,
                )
                actions.append(data.data)
                actions[-1]["derived"] = data.derived
            self.save_bulk(actions)

    def upload(
        self,
        file: Path,
        extra_infos: Union[Iterable, dict, None] = None,
        store_calc: bool = True,
    ):
        """
        Upload data from a file to the database.

        Parameters
        ----------
        file: Path
            Path to file to be uploaded
        extra_infos: Union[Iterable, dict, None], optional
            Extra information to store in the document with the atoms data.
            Default is `None`.
        store_calc: bool, optional
            Whether to store data from the calculator attached to atoms.
            Default is `True`.
        """

        if isinstance(file, str):
            file = Path(file)

        extra_info = {}
        if extra_infos:
            for info in extra_infos:
                extra_info.update(extras.parser.parse(info))  # type: ignore

        extra_info["filename"] = str(file)

        data = iread(str(file))
        self.push(data, extra_info, store_calc=store_calc)

    def get_items(
        self, query: Union[dict, str, None] = None
    ) -> Generator[dict, None, None]:
        """
        Get data as a dictionary from documents in the database.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to get data from. Default is `None`.

        Returns
        -------
        Generator for dictionary of data.
        """
        query = self.parser(query)
        query = {
            "query": query,
        }

        for hit in helpers.scan(
            self.client,
            index=self.index_name,
            query=query,
        ):
            yield {"_id": hit["_id"], **hit["_source"]}

    def get_atoms(
        self, query: Union[dict, str, None] = None
    ) -> Generator[Atoms, None, None]:
        """
        Get data as Atoms object from documents in the database.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to get data from. Default is `None`.

        Returns
        -------
        Generator for AtomsModel object of data.
        """
        query = self.parser(query)
        query = {
            "query": query,
        }

        for hit in helpers.scan(
            self.client,
            index=self.index_name,
            query=query,
        ):
            yield AtomsModel(None, None, hit["_source"]).to_ase()

    def count(self, query: Union[dict, str, None] = None) -> int:
        """
        Counts number of documents in the database.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to be counted. Default is `None`.

        Returns
        -------
        Count of number of documents.
        """
        logger.info("query; {}".format(query))
        query = self.parser(query)
        body = {
            "query": query,
        }

        return self.client.count(index=self.index_name, body=body)["count"]

    def property(self, name, query: Union[dict, str, None] = None) -> list:
        """
        Gets all values of a specified property for matching documents in the
        database. This method is very slow, so it is preferable to use
        alternative methods where possible, such as count_property.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to get properties from. Default is `None`.

        Returns
        -------
        List of values for the specified property for all matching documents.
        """
        query = self.parser(query)
        query = {
            "query": query,
        }

        return [
            hit["_source"][format(name)]
            for hit in helpers.scan(
                self.client,
                index=self.index_name,
                query=query,
                stored_fields=format(name),
                _source=format(name),
            )
        ]

    def count_property(self, name, query: Union[dict, str, None] = None) -> dict:
        """
        Counts values of a specified property for matching documents in the
        database. This method much faster than performing a Count on the list
        returned by self.property, so this method should be used preferentially.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to count properties from. Default is `None`.

        Returns
        -------
        Dictionary of values and counts for the specified property for all
        matching documents.
        """
        query = self.parser(query)

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

        for val in self.client.search(index=self.index_name, body=body,)[
            "aggregations"
        ][format(name)]["buckets"]:
            prop[val["key"]] = val["doc_count"]

        return prop

    def properties(self, query: Union[dict, str, None] = None) -> dict:
        """
        Gets lists of all properties from matching documents, separated into
        info, derived, and array properties.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to get properties from. Default is `None`.

        Returns
        -------
        Dictionary of properties, with keys corresponding to info, derived,
        and arrays of properties, and values corresponding to a list of
        the properties of that type.
        """
        query = self.parser(query)

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

            derived = ["info_keys", "derived_keys", "arrays_keys"]
            for label in derived:
                count = res["aggregations"][label]["doc_count"]
                if count > 0:
                    key = label.split("_")[0]
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
        # TODO: Probably it would be nicer to store the type info in the database from the beginning.
        atoms = self.client.search(
            index=self.index_name,
            body={"size": 1, "query": {"exists": {"field": prop}}},
        )

        data = atoms["hits"]["hits"][0]["_source"][prop]

        if category == "arrays":
            if type(data[0]) == list:
                return "array({}, N x {})".format(
                    map_types[type(data[0][0])], len(data[0])
                )
            else:
                return "vector({}, N)".format(map_types[type(data[0])])

        if type(data) == list:
            if type(data[0]) == list:
                if type(data[0][0]) == list:
                    return "list(list(...)"
                else:
                    return "array({})".format(map_types[type(data[0][0])])
            else:
                return "vector({})".format(map_types[type(data[0])])
        else:
            return "scalar({})".format(map_types[type(data)])

    def count_properties(self, query: Union[dict, str, None] = None) -> dict:
        """
        Counts all properties from matching documents.

        Parameters
        ----------
        query: Union[dict, str, None]
            Query to filter documents to count properties from. Default is `None`.

        Returns
        -------
        Dictionary of properties, with keys property names, and values
        corresponding to their counts, categories and data types.
        """
        query = self.parser(query)
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

            derived = ["info_keys", "derived_keys", "arrays_keys"]
            for label in derived:
                count = res["aggregations"][label]["doc_count"]
                if count > 0:
                    properties[key] = {
                        "count": count,
                        "category": label.split("_")[0],
                        "dtype": self.get_type_of_property(key, label.split("_")[0]),
                    }

        return properties

    def add_property(self, data: dict, query: Union[dict, str, None] = None):
        """
        Adds properties to matching documents.

        Parameters
        ----------
        data: dict
            Property key-value pairs to be added to matching documents.
        query: Union[dict, str, None]
            Query to filter documents to add properties to. Default is `None`.
        """
        logger.info("add: data={}, query={}".format(data, query))
        query = self.parser(query)

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
        self, name: str, new_name: str, query: Union[dict, str, None] = None
    ):
        """
        Renames property for all matching documents.

        Parameters
        ----------
        name: str
            Current name of property to be renamed.
        new_name: str
            New name of property to be renamed.
        query: Union[dict, str, None]
            Query to filter documents to rename property. Default is `None`.
        """
        logger.info("rename: query={}, old={}, new={}".format(query, name, new_name))
        query = self.parser(query)

        script_txt = f"if (!ctx._source.containsKey('{new_name}')) {{ "
        script_txt += (
            f"ctx._source.{new_name} = ctx._source.{name};"
            " ctx._source.remove('params.name');"
        )

        script_txt += f"for (int i=0; i<ctx._source.derived.info_keys.length; i++) {{"
        script_txt += f"if (ctx._source.derived.info_keys[i] == params.name) {{ "
        script_txt += f"ctx._source.derived.info_keys[i] = params.new_name;}}}}}}"

        body = {
            "script": {
                "source": script_txt,
                "lang": "painless",
                "params": {"name": name, "new_name": new_name},
            },
            "query": query,
        }

        self.client.update_by_query(index=self.index_name, body=body)

    def delete_property(self, name: str, query: Union[dict, str, None] = None):
        """
        Deletes property from all matching documents.

        Parameters
        ----------
        name: str
            Name of property to be deleted from documents.
        query: Union[dict, str, None]
            Query to filter documents to have property deleted. Default is `None`.
        """
        logger.info("delete: query={}, porperty={}".format(name, query))
        query = self.parser(query)

        script_txt = f"if (ctx._source.containsKey('{name}')) {{ "
        script_txt += f"ctx._source.remove('params.name');"
        script_txt += f"for (int i=0; i<ctx._source.derived.info_keys.length; i++) {{"
        script_txt += f"if (ctx._source.derived.info_keys[i] == params.name) {{ "
        script_txt += f"ctx._source.derived.info_keys.remove(i);}}}}}}"

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
        self, name: str, query: Union[dict, str, None] = None, **kwargs
    ) -> Union[dict, None]:
        """
        Calculate histogram statistics for a property from all matching documents.

        Parameters
        ----------
        name: str
            Name of property.
        query: Union[dict, str, None]
            Query to filter documents. Default is `None`.

        Returns
        -------
        Dictionary containing histogram statistics, including the number of
        bins, edges, counts, min, max, and standard deviation.
        """
        query = self.parser(query)

        data = self.property(name, query)
        return histogram(name, data, **kwargs)

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
            "{}(".format(self.__class__.__name__)
            + "url={}:{}, ".format(host, port)
            + "index={}) ".format(self.index_name)
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


def histogram(name, data, **kwargs):
    if not data:
        return None

    elif data and isinstance(data, list):

        ptype = type(data[0])

        if not all(isinstance(x, ptype) for x in data):
            print("Mixed type error of the {} property!".format(name))
            return None

        if ptype == float:
            bins = kwargs.get("bins", 10)
            return _hist_float(name, data, bins)

        elif ptype == int:
            bins = kwargs.get("bins", 10)
            return _hist_int(name, data, bins)

        elif ptype == str:
            return _hist_str(name, data, **kwargs)

        elif ptype == datetime:
            bins = kwargs.get("bins", 10)
            return _hist_date(name, data, bins)

        else:
            print(
                "{}: Histogram for list of {} types are not supported!".format(
                    name, type(data[0])
                )
            )
            logger.info(
                "{}: Histogram for list of {} types are not supported!".format(
                    name, type(data[0])
                )
            )

    else:
        logger.info(
            "{}: Histogram for {} types are not supported!".format(name, type(data))
        )
        return None


def _hist_float(name, data, bins=10):
    data = np.array(data)
    hist, bin_edges = np.histogram(data, bins=bins)

    return {
        "type": "hist_float",
        "name": name,
        "bins": bins,
        "edges": bin_edges,
        "counts": hist,
        "min": data.min(),
        "max": data.max(),
        "median": data.mean(),
        "std": data.std(),
        "var": data.var(),
    }


def _hist_date(name, data, bins=10):
    hist_data = np.array([t.timestamp() for t in data])
    hist, bin_edges = np.histogram(hist_data, bins=bins)

    fromtimestamp = datetime.fromtimestamp

    return {
        "type": "hist_date",
        "name": name,
        "bins": bins,
        "edges": [fromtimestamp(d) for d in bin_edges],
        "counts": hist,
        "min": fromtimestamp(hist_data.min()),
        "max": fromtimestamp(hist_data.max()),
        "median": fromtimestamp(hist_data.mean()),
        "std": fromtimestamp(hist_data.std()),
        "var": fromtimestamp(hist_data.var()),
    }


def _hist_int(name, data, bins=10):
    data = np.array(data)
    delta = max(data) - min(data) + 1

    if bins > delta:
        bins = delta

    hist, bin_edges = np.histogram(data, bins=bins)

    return {
        "type": "hist_int",
        "name": name,
        "bins": bins,
        "edges": bin_edges,
        "counts": hist,
        "min": data.min(),
        "max": data.max(),
        "median": data.mean(),
        "std": data.std(),
        "var": data.var(),
    }


def _hist_str(name, data, bins=10, truncate=20):
    n_unique = len(set(data))

    if truncate:
        # data = (item[:truncate] for item in data)
        data = (
            item[:truncate] + "..." if len(item) > truncate else item for item in data
        )

    data = Counter(data)

    if bins:
        labels, counts = zip(*sorted(data.items(), key=itemgetter(1, 0), reverse=True))
    else:
        labels, counts = zip(*data.items())

    return {
        "type": "hist_str",
        "name": name,
        "total": sum(data.values()),
        "unique": n_unique,
        "labels": labels[:bins],
        "counts": counts[:bins],
    }


if __name__ == "__main__":
    db = OpenSearchDatabase(username="admin", password="admin")
    print(db.info())
