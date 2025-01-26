# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar
# This program is distributed under the MIT License, see LICENSE.md.

import logging

from abcd import ABCD

from abcd.frontends.commandline.config import Config
from abcd.parsers.queries import parser

logger = logging.getLogger(__name__)


def init_config(func):
    config = Config.load()

    def wrapper(*args, **kwargs):
        func(*args, config=config, **kwargs)

    return wrapper


def init_db(func):
    def wrapper(*args, config, **kwargs):
        url = config.get("url", None)

        if url is None:
            print("Please use abcd login first!")
            exit(1)

        db = ABCD.from_url(url=url)

        # TODO: AST.from_string() ?!
        # TODO: parser should accept list
        # TODO: better ast optimisation

        query_list = []
        for q in kwargs.pop("default_query", []):
            query_list.append(parser(q))

        for q in kwargs.pop("query", []):
            query_list.append(parser(q))

        if not query_list:
            query = None
        elif len(query_list) == 1:
            query = query_list[0]
        else:
            query = ("AND", *query_list)

        func(*args, db=db, query=query, **kwargs)

    return wrapper


def check_remote(func):
    def wrapper(*args, **kwargs):
        if kwargs.pop("remote"):
            print("In read only mode, you can't modify the data in the database")
            exit(1)

        func(*args, **kwargs)

    return wrapper
