import logging
import functools
from abcd import ABCD
from abcd.frontends.commandline.config import Config

logger = logging.getLogger(__name__)


def init_config(func):
    config = Config.load()

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        func(*args, config=config, **kwargs)

    return wrapper


def init_db(func):
    @functools.wraps(func)
    def wrapper(*args, config, **kwargs):
        url = config.get("url", None)
        use_ssl = config.get("use_ssl", None)

        if url is None:
            print("Please use abcd login first!")
            exit(1)

        if use_ssl is None:
            print("use_ssl has not been saved. Please login again")
            exit(1)

        db = ABCD.from_url(url=url, use_ssl=use_ssl)

        # TODO: AST.from_string() ?!
        # TODO: better ast optimisation

        query_list = []
        for q in kwargs.pop("default_query", []):
            query_list.append(q)

        for q in kwargs.pop("query", []):
            query_list.append(q)

        func(*args, db=db, query=query_list, **kwargs)

    return wrapper


def check_remote(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if kwargs.pop("remote"):
            print("In read only mode, you can't modify the data in the database")
            exit(1)

        func(*args, **kwargs)

    return wrapper
