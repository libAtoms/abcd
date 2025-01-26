# Copyright (c) 2025.
# Authors: Ádám Fekete
# This program is distributed under the MIT License, see LICENSE.md.

import logging
from abc import ABCMeta

logger = logging.getLogger(__name__)


class Query(dict):
    pass


class AbstractQuerySet(metaclass=ABCMeta):
    pass
