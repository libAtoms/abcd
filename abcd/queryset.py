# Copyright (c) 2025.
# Authors: Ádám Fekete
# This program is distributed under the MIT License, see LICENSE.md.

from abc import ABCMeta
import logging

logger = logging.getLogger(__name__)


class Query(dict):
    pass


class AbstractQuerySet(metaclass=ABCMeta):  # noqa: B024
    pass
