import types
import logging
import numpy as np

from typing import Union, Iterable
from os import linesep
from operator import itemgetter
from collections import Counter
from datetime import datetime

from ase import Atoms
from ase.io import iread

import abcd.errors
from abcd.model import AbstractModel
from abcd.database import AbstractABCD
from abcd.queryset import AbstractQuerySet
from abcd.parsers import extras

import pymongo.errors
from pymongo import MongoClient
from bson import ObjectId

from pathlib import Path

logger = logging.getLogger(__name__)

map_types = {
    bool: "bool",
    float: "float",
    int: "int",
    str: "str",
    datetime: "date",
    dict: "dict",
}


class AtomsModel(AbstractModel):
    def __init__(self, collection=None, dict=None):
        super().__init__(dict)

        self._collection = collection

    @classmethod
    def from_atoms(cls, collection, atoms: Atoms, extra_info=None, store_calc=True):
        obj = super().from_atoms(atoms, extra_info, store_calc)
        obj._collection = collection
        return obj

    @property
    def _id(self):
        return self.get("_id", None)

    def save(self):
        if not self._id:
            self._collection.insert_one(self)
        else:
            new_values = {"$set": self}
            self._collection.update_one({"_id": ObjectId(self._id)}, new_values)

    def remove(self):
        if self._id:
            self._collection.remove({"_id": ObjectId(self._id)})
            self.clear()


class MongoQuery(AbstractQuerySet):
    def __init__(self):
        pass

    def visit(self, syntax_tree):
        op, *args = syntax_tree
        try:
            fun = self.__getattribute__("visit_" + op.lower())
            return fun(*args)
        except KeyError:
            pass

    def visit_name(self, field):
        return {field: {"$exists": True}}

    def visit_not(self, value):
        _, field = value
        return {field: {"$exists": False}}

    def visit_and(self, *args):
        print(args)
        return {"$and": [self.visit(arg) for arg in args]}
        # TODO recursively combining all the and statements
        # out = {}
        # for arg in args:
        #     a = self.visit(arg)
        #
        #     out.update(**a)
        # return out

    def visit_or(self, *args):
        return {"$or": [self.visit(arg) for arg in args]}

    def visit_eq(self, field, value):
        return {field[1]: value[1]}

    def visit_re(self, field, value):
        return {field[1]: {"$regex": value[1]}}

    def visit_gt(self, field, value):
        return {field[1]: {"$gt": value[1]}}

    def visit_gte(self, field, value):
        return {field[1]: {"$gte": value[1]}}

    def visit_lt(self, field, value):
        return {field[1]: {"$lt": value[1]}}

    def visit_lte(self, field, value):
        return {field[1]: {"$lte": value[1]}}

    def visit_in(self, field, *values):
        return {field[1]: {"$in": [value[1] for value in values]}}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __call__(self, ast):
        logger.info("parsed ast: {}".format(ast))

        if isinstance(ast, dict):
            return ast
        elif isinstance(ast, str):
            from abcd.parsers.queries import parser

            p = parser(ast)
            return self.visit(p)

        elif isinstance(ast, list):
            from abcd.parsers.queries import parser

            if len(ast) == 0:
                return {}
            else:
                ast = ("AND", *[parser(q) for q in ast])
                return self.visit(ast)

        return self.visit(ast) if ast else {}


parser = MongoQuery()


def parse_query(func):
    def wrapper(*args, query=None, **kwargs):
        print(func)
        print((args, query, kwargs))
        query = parser(query)

        func(*args, **kwargs, query=query)

    return wrapper


class MongoDatabase(AbstractABCD):
    """Wrapper to make database operations easy"""

    def __init__(
        self,
        host="localhost",
        port=27017,
        db_name="abcd",
        collection_name="atoms",
        username=None,
        password=None,
        authSource="admin",
        **kwargs
    ):
        super().__init__()

        logger.info(
            (
                host,
                port,
                db_name,
                collection_name,
                username,
                password,
                authSource,
                kwargs,
            )
        )

        self.client = MongoClient(
            host=host,
            port=port,
            username=username,
            password=password,
            authSource=authSource,
        )

        try:
            info = self.client.server_info()  # Forces a call.
            logger.info("DB info: {}".format(info))

        except pymongo.errors.OperationFailure:
            raise abcd.errors.AuthenticationError()

        except pymongo.errors.ServerSelectionTimeoutError:
            raise abcd.errors.TimeoutError()

        self.db = self.client[db_name]
        self.collection = self.db[collection_name]

    def info(self):
        host, port = self.client.address

        return {
            "host": host,
            "port": port,
            "db": self.db.name,
            "collection": self.collection.name,
            "number of confs": self.collection.count_documents({}),
            "type": "mongodb",
        }

    def delete(self, query=None):
        query = parser(query)
        return self.collection.delete_many(query)

    def destroy(self):
        self.collection.drop()

    def push(self, atoms: Union[Atoms, Iterable], extra_info=None, store_calc=True):
        if extra_info and isinstance(extra_info, str):
            extra_info = extras.parser.parse(extra_info)

        if isinstance(atoms, Atoms):
            data = AtomsModel.from_atoms(
                self.collection, atoms, extra_info=extra_info, store_calc=store_calc
            )
            data.save()
            # self.collection.insert_one(data)

        elif isinstance(atoms, types.GeneratorType) or isinstance(atoms, list):
            for i, item in enumerate(atoms):
                if isinstance(extra_info, list):
                    info = extra_info[i]
                else:
                    info = extra_info
                data = AtomsModel.from_atoms(
                    self.collection, item, extra_info=info, store_calc=store_calc
                )
                data.save()

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

    def get_items(self, query=None):
        # TODO: better method for aggregations
        query = parser(query)
        for dct in self.db.atoms.find(query):
            yield dct

    def get_atoms(self, query=None):
        query = parser(query)
        for dct in self.db.atoms.find(query):
            yield AtomsModel(None, dct).to_ase()

    def count(self, query=None):
        query = parser(query)
        logger.info("query; {}".format(query))

        if not query:
            query = {}

        return self.db.atoms.count_documents(query)

    def property(self, name, query=None):
        query = parser(query)

        pipeline = [
            {"$match": query},
            {"$match": {"{}".format(name): {"$exists": True}}},
            {"$project": {"_id": False, "data": "${}".format(name)}},
        ]

        return [val["data"] for val in self.db.atoms.aggregate(pipeline)]

    def properties(self, query=None):
        query = parser(query)
        properties = {}

        pipeline = [
            {"$match": query},
            {"$unwind": "$derived.info_keys"},
            {"$group": {"_id": "$derived.info_keys"}},
        ]
        properties["info"] = [
            value["_id"] for value in self.db.atoms.aggregate(pipeline)
        ]

        pipeline = [
            {"$match": query},
            {"$unwind": "$derived.arrays_keys"},
            {"$group": {"_id": "$derived.arrays_keys"}},
        ]
        properties["arrays"] = [
            value["_id"] for value in self.db.atoms.aggregate(pipeline)
        ]

        pipeline = [
            {"$match": query},
            {"$unwind": "$derived.derived_keys"},
            {"$group": {"_id": "$derived.derived_keys"}},
        ]
        properties["derived"] = [
            value["_id"] for value in self.db.atoms.aggregate(pipeline)
        ]

        return properties

    def get_type_of_property(self, prop, category):
        # TODO: Probably it would be nicer to store the type info in the database from the beginning.
        atoms = self.db.atoms.find_one({prop: {"$exists": True}})
        data = atoms[prop]

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

    def count_properties(self, query=None):
        query = parser(query)

        properties = {}

        pipeline = [
            {"$match": query},
            {"$unwind": "$derived.info_keys"},
            {"$group": {"_id": "$derived.info_keys", "count": {"$sum": 1}}},
        ]

        info_keys = self.db.atoms.aggregate(pipeline)
        for val in info_keys:
            properties[val["_id"]] = {
                "count": val["count"],
                "category": "info",
                "dtype": self.get_type_of_property(val["_id"], "info"),
            }

        pipeline = [
            {"$match": query},
            {"$unwind": "$derived.arrays_keys"},
            {"$group": {"_id": "$derived.arrays_keys", "count": {"$sum": 1}}},
        ]
        arrays_keys = list(self.db.atoms.aggregate(pipeline))
        for val in arrays_keys:
            properties[val["_id"]] = {
                "count": val["count"],
                "category": "arrays",
                "dtype": self.get_type_of_property(val["_id"], "arrays"),
            }

        pipeline = [
            {"$match": query},
            {"$unwind": "$derived.derived_keys"},
            {"$group": {"_id": "$derived.derived_keys", "count": {"$sum": 1}}},
        ]
        arrays_keys = list(self.db.atoms.aggregate(pipeline))
        for val in arrays_keys:
            properties[val["_id"]] = {
                "count": val["count"],
                "category": "derived",
                "dtype": self.get_type_of_property(val["_id"], "derived"),
            }

        return properties

    def add_property(self, data, query=None):
        logger.info("add: data={}, query={}".format(data, query))

        self.collection.update_many(
            parser(query),
            {
                "$push": {"derived.info_keys": {"$each": list(data.keys())}},
                "$set": data,
            },
        )

    def rename_property(self, name, new_name, query=None):
        logger.info("rename: query={}, old={}, new={}".format(query, name, new_name))
        # TODO name in derived.info_keys OR name in derived.arrays_keys OR name in derived.derived_keys
        self.collection.update_many(
            parser(query), {"$push": {"derived.info_keys": new_name}}
        )

        self.collection.update_many(
            parser(query),
            {"$pull": {"derived.info_keys": name}, "$rename": {name: new_name}},
        )

        # self.collection.update_many(
        #     parser(query + ['arrays.{}'.format(name)]),
        #     {'$push': {'derived.arrays_keys': new_name}})
        #
        # self.collection.update_many(
        #     parser(query + ['arrays.{}'.format(name)]),
        #     {'$pull': {'derived.arrays_keys': name},
        #      '$rename': {'arrays.{}'.format(name): 'arrays.{}'.format(new_name)}})

    def delete_property(self, name, query=None):
        logger.info("delete: query={}, porperty={}".format(name, query))

        self.collection.update_many(
            parser(query),
            {
                "$pull": {"derived.info_keys": name, "derived.arrays_keys": name},
                "$unset": {name: ""},
            },
        )

    def hist(self, name, query=None, **kwargs):
        data = self.property(name, query)
        return histogram(name, data, **kwargs)

    def exec(self, code, query=None):
        # TODO: Separate python environment with its own packages loaded

        for dct in self.get_items(query):
            atoms = AtomsModel(self.collection, dct)
            exec(code)

    def __repr__(self):
        host, port = self.client.address

        return (
            "{}(".format(self.__class__.__name__)
            + "url={}:{}, ".format(host, port)
            + "db={}, ".format(self.db.name)
            + "collection={})".format(self.collection.name)
        )

    def _repr_html_(self):
        """Jupyter notebook representation"""
        return "<b>ABCD MongoDB database</b>"

    def print_info(self):
        """shows basic information about the connected database"""

        out = linesep.join(
            [
                "{:=^50}".format(" ABCD MongoDB "),
                "{:>10}: {}".format("type", "mongodb"),
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
    # import json
    # from ase.io import iread
    # from pprint import pprint
    # from server.styles.myjson import JSONEncoderOld, JSONDecoderOld, JSONEncoder

    print("hello")
    db = MongoDatabase(username="mongoadmin", password="secret")
    print(db.info())
    print(db.count())

    print(db.hist("uploaded"))

    # for atoms in iread('../../tutorials/data/bcc_bulk_54_expanded_2_high.xyz', index=slice(None)):
    #     # print(at)
    #     atoms.calc.results['forces'] = atoms.arrays['force']
    #     # at.arrays['force'] = None
    #
    #     json_data = json.dumps(atoms, cls=JSONEncoderOld)
    #     print(json_data)
    #
    #     atom_dict = json.loads(json_data, cls=JSONDecoderOld)
    #     print(atom_dict)
    #
    #     print(atoms == atom_dict)
    #
    # with JSONEncoder() as encoder:
    #     data = encoder.encode(atoms)
    #
    # print(data)
    #
    # with DictEncoder() as encoder:
    #     data = encoder.encode(atoms)
    #
    # pprint(data)
    #
    # db.collection.insert_one(DictEncoder().encode(atoms))
