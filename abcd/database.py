# Copyright (c) 2025.
# Authors: Ádám Fekete
# This program is distributed under the MIT License, see LICENSE.md.

# ruff: noqa: B024, B027
from abc import ABCMeta
import logging

logger = logging.getLogger(__name__)


class AbstractABCD(metaclass=ABCMeta):
    """Factory method"""

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
