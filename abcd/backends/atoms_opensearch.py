import types
import logging

from typing import Union, Iterable
from os import linesep
from datetime import datetime
from collections import Counter
from operator import itemgetter


import numpy as np

from ase import Atoms
from ase.io import iread

import abcd.errors
from abcd.model import AbstractModel
from abcd.database import AbstractABCD
from abcd.queryset import AbstractQuerySet
from abcd.parsers import extras

from pathlib import Path

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
    dict: "dict"
}


class OpenSearchQuery(AbstractQuerySet):

    def __init__(self, client, index_name, analyse_schema=True):
        if analyse_schema:
            schema = client.indices.get_mapping()[index_name]
            schema_analizer = SchemaAnalyzer(schema)
            self.query_builder = ElasticsearchQueryBuilder(**schema_analizer.query_builder_options())
        else:
            self.query_builder = ElasticsearchQueryBuilder()

    def __call__(self, ast):
        logger.info('parsed ast: {}'.format(ast))

        if isinstance(ast, dict):
            return ast
        elif isinstance(ast, str):
            tree = parser.parse(ast)
            return self.query_builder(tree)

        return ast if ast else None


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

    def save(self):
        if not self._id:
            body = {}
            body.update(self.data)
            body["derived"] = self.derived
            self._client.index(index=self._index_name, body=body)


class OpenSearchDatabase(AbstractABCD):
    """Wrapper to make database operations easy"""

    def __init__(
            self,
            host="localhost",
            port=9200,
            db="abcd",
            index_name="atoms",
            username="admin",
            password="admin",
            analyse_schema=True,
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

        self.db = db
        self.index_name = index_name
        self.create()
        self.parser = OpenSearchQuery(self.client, self.index_name, analyse_schema)

    def info(self):
        host = self.client.transport.hosts[0]["host"]
        port = self.client.transport.hosts[0]["port"]

        self.client.indices.refresh(index=self.index_name)
        return {
            "host": host,
            "port": port,
            "db": self.db,
            "index": self.index_name,
            "number of confs": self.client.count(index=self.index_name)["count"],
            "type": "opensearch"
        }

    def delete(self, query=None):
        query = self.parser(query)
        if not query:
            query = {
                "match_all": {}
            }

        self.client.delete_by_query(
            index=self.index_name,
            body={
                "query": query,
            },
        )

    def destroy(self):
        self.client.indices.delete(index=self.index_name, ignore=404)

    def create(self):
        self.client.indices.create(index=self.index_name, ignore=400)

    def save_bulk(self, actions: Iterable):
        helpers.bulk(client=self.client, actions=actions, index=self.index_name)

    def push(self, atoms: Union[Atoms, Iterable], extra_info=None, store_calc=True):

        if extra_info and isinstance(extra_info, str):
            extra_info = extras.parser.parse(extra_info)

        # Could combine into single data.save, but keep separate for option of bulk insertion?
        if isinstance(atoms, Atoms):
            data = AtomsModel.from_atoms(self.client, self.index_name, atoms, extra_info=extra_info, store_calc=store_calc)
            data.save()

        elif isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):

            actions = []
            for item in atoms:
                data = AtomsModel.from_atoms(self.client, self.index_name, item, extra_info=extra_info, store_calc=store_calc)
                actions.append(data.data)
                actions[-1]["derived"] = data.derived
            self.save_bulk(actions)

    def upload(self, file: Path, extra_infos=None, store_calc=True):

        if isinstance(file, str):
            file = Path(file)

        extra_info = {}
        if extra_infos:
            for info in extra_infos:
                extra_info.update(extras.parser.parse(info))

        extra_info["filename"] = str(file)

        data = iread(str(file))
        self.push(data, extra_info, store_calc=store_calc)

    def get_atoms(self, query=None):
        query = self.parser(query)
        if not query:
            query = {
                "query": {
                    "match_all": {}
                }
            }

        for hit in helpers.scan(
            self.client,
            index=self.index_name,
            query=query,
        ):
            yield AtomsModel(None, None, hit["_source"]).to_ase()

    def count(self, query=None):
        query = self.parser(query)
        logger.info("query; {}".format(query))

        if not query:
            query = {
                "match_all": {}
            }

        return self.client.count(index=self.index_name, body={"query": query})["count"]

    # Slow - use count_property where possible!
    def property(self, name, query=None):
        query = self.parser(query)
        if not query:
            query = {
                "match_all": {}
            }

        body = {
            "query": query,
        }

        return [hit["_source"][format(name)] for hit in helpers.scan(
            self.client,
            index=self.index_name,
            query=body,
            stored_fields=format(name),
            _source=format(name),
        )]

    def count_property(self, name, query=None):
        query = self.parser(query)
        if not query:
            query = {
                "match_all": {}
            }

        body = {
            "size" : 0,
            "query": query,
            "aggs": {
                format(name): {
                    "terms": {
                        "field": format(name),
                        "size": 10000, # Use composite for all results?
                    },
                },
            }
        }

        prop = {}

        for val in self.client.search(
            index=self.index_name,
            body=body,
        )["aggregations"][format(name)]["buckets"]:
            prop[val["key"]] = val["doc_count"]

        return prop

    def properties(self, query=None):
        query = self.parser(query)
        if not query:
            query = {
                "match_all": {}
            }

        properties = {}

        for prop in self.client.indices.get_mapping(self.index_name)[self.index_name]["mappings"]["properties"].keys():

            body = {
                "size" : 0,
                "query": query,
                "aggs": {
                    "info_keys": {
                        "filter": {
                            "term": {
                                "derived.info_keys.keyword": prop
                            }
                        },
                    },
                    "derived_keys": {
                        "filter": {
                            "term": {
                                "derived.derived_keys.keyword": prop
                            }
                        },
                    },
                    "arrays_keys": {
                        "filter": {
                            "term": {
                                "derived.arrays_keys.keyword": prop
                            }
                        },
                    },
                }
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

    def get_type_of_property(self, prop, category):
        # TODO: Probably it would be nicer to store the type info in the database from the beginning.
        atoms = self.client.search(
            index=self.index_name,
            body = {
                "size" : 1,
                "query": {
                   "exists" : {
                        "field": prop
                    }
                }
            }
        )

        data = atoms["hits"]["hits"][0]["_source"][prop]

        if category == "arrays":
            if type(data[0]) == list:
                return "array({}, N x {})".format(map_types[type(data[0][0])], len(data[0]))
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

    def count_properties(self, query=None):
        query = self.parser(query)
        if not query:
            query = {
                "match_all": {}
            }

        properties = {}

        try:
            keys = self.client.indices.get_mapping(self.index_name)[self.index_name]["mappings"]["properties"].keys()
        except KeyError:
            return properties

        for key in keys:

            body = {
                "size" : 0,
                "query": query,
                "aggs": {
                    "info_keys": {
                        "filter": {
                            "term": {
                                "derived.info_keys.keyword": key
                            }
                        },
                    },
                    "derived_keys": {
                        "filter": {
                            "term": {
                                "derived.derived_keys.keyword": key
                            }
                        },
                    },
                    "arrays_keys": {
                        "filter": {
                            "term": {
                                "derived.arrays_keys.keyword": key
                            }
                        },
                    },
                }
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
                        "dtype": self.get_type_of_property(key, label.split("_")[0])
                    }

        return properties

    def add_property(self, data, query=None):
        logger.info('add: data={}, query={}'.format(data, query))
        query = self.parser(query)

        script_txt = "ctx._source.derived.info_keys.addAll(params.keys);"
        for key, val in data.items():
            script_txt += f"ctx._source.{key} = '{val}';"

        body = {
            "script": {
                "source": script_txt,
                "lang": "painless",
                "params" : {
                    "keys" : list(data.keys())
                },
            },
            "query": query
        }

        self.client.update_by_query(
            index=self.index_name,
            body=body,
        )

    def rename_property(self, name, new_name, query=None):
        logger.info('rename: query={}, old={}, new={}'.format(query, name, new_name))
        query = self.parser(query)

        script_txt = f"if (!ctx._source.containsKey('{new_name}')) {{ "
        script_txt += f"ctx._source.{new_name} = ctx._source.{name}; ctx._source.remove('params.name');"

        script_txt += f"for (int i=0; i<ctx._source.derived.info_keys.length; i++) {{"
        script_txt += f"if (ctx._source.derived.info_keys[i] == params.name) {{ "
        script_txt += f"ctx._source.derived.info_keys[i] = params.new_name;}}}}}}"

        body = {
            "script": {
                "source": script_txt,
                "lang": "painless",
                "params": {
                    "name": name,
                    "new_name": new_name
                },
            },
            "query": query
        }

        self.client.update_by_query(
            index=self.index_name,
            body=body
        )

    def delete_property(self, name, query=None):
        logger.info('delete: query={}, porperty={}'.format(name, query))
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
            "query": query
        }

        self.client.update_by_query(
            index=self.index_name,
            body=body
        )

    def hist(self, name, query=None, **kwargs):
        query = self.parser(query)

        data = self.property(name, query)
        return histogram(name, data, **kwargs)

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


def histogram(name, data, **kwargs):
    if not data:
        return None

    elif data and isinstance(data, list):

        ptype = type(data[0])

        if not all(isinstance(x, ptype) for x in data):
            print("Mixed type error of the {} property!".format(name))
            return None

        if ptype == float:
            bins = kwargs.get('bins', 10)
            return _hist_float(name, data, bins)

        elif ptype == int:
            bins = kwargs.get('bins', 10)
            return _hist_int(name, data, bins)

        elif ptype == str:
            return _hist_str(name, data, **kwargs)

        elif ptype == datetime:
            bins = kwargs.get('bins', 10)
            return _hist_date(name, data, bins)

        else:
            print('{}: Histogram for list of {} types are not supported!'.format(name, type(data[0])))
            logger.info('{}: Histogram for list of {} types are not supported!'.format(name, type(data[0])))

    else:
        logger.info('{}: Histogram for {} types are not supported!'.format(name, type(data)))
        return None


def _hist_float(name, data, bins=10):
    data = np.array(data)
    hist, bin_edges = np.histogram(data, bins=bins)

    return {
        'type': 'hist_float',
        'name': name,
        'bins': bins,
        'edges': bin_edges,
        'counts': hist,
        'min': data.min(),
        'max': data.max(),
        'median': data.mean(),
        'std': data.std(),
        'var': data.var()
    }


def _hist_date(name, data, bins=10):
    hist_data = np.array([t.timestamp() for t in data])
    hist, bin_edges = np.histogram(hist_data, bins=bins)

    fromtimestamp = datetime.fromtimestamp

    return {
        'type': 'hist_date',
        'name': name,
        'bins': bins,
        'edges': [fromtimestamp(d) for d in bin_edges],
        'counts': hist,
        'min': fromtimestamp(hist_data.min()),
        'max': fromtimestamp(hist_data.max()),
        'median': fromtimestamp(hist_data.mean()),
        'std': fromtimestamp(hist_data.std()),
        'var': fromtimestamp(hist_data.var())
    }


def _hist_int(name, data, bins=10):
    data = np.array(data)
    delta = max(data) - min(data) + 1

    if bins > delta:
        bins = delta

    hist, bin_edges = np.histogram(data, bins=bins)

    return {
        'type': 'hist_int',
        'name': name,
        'bins': bins,
        'edges': bin_edges,
        'counts': hist,
        'min': data.min(),
        'max': data.max(),
        'median': data.mean(),
        'std': data.std(),
        'var': data.var()
    }


def _hist_str(name, data, bins=10, truncate=20):
    n_unique = len(set(data))

    if truncate:
        # data = (item[:truncate] for item in data)
        data = (item[:truncate] + '...' if len(item) > truncate else item for item in data)

    data = Counter(data)

    if bins:
        labels, counts = zip(*sorted(data.items(), key=itemgetter(1, 0), reverse=True))
    else:
        labels, counts = zip(*data.items())

    return {
        'type': 'hist_str',
        'name': name,
        'total': sum(data.values()),
        'unique': n_unique,
        'labels': labels[:bins],
        'counts': counts[:bins]
    }


if __name__ == "__main__":
    db = OpenSearchDatabase(username="admin", password="admin")
    print(db.info())
