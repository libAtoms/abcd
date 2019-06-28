from abc import ABCMeta, abstractmethod
from enum import Enum


class ConnectionType(Enum):
    mongodb = 1
    http = 2


class URLError(Exception):
    pass


class AuthenticationError(Exception):
    pass


class Database(object, metaclass=ABCMeta):

    @abstractmethod
    def __init__(self):
        pass

    def info(self):
        pass

    def push(self, atoms):
        pass

    def pull(self, query=None, properties=None):
        pass

    def query(self, query_string):
        pass

    def destroy(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __repr__(self):
        pass

    def _repr_html_(self):
        pass

    def print_info(self):
        pass
