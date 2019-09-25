import logging
from abc import ABCMeta

logger = logging.getLogger(__name__)


class Query(dict):
    pass


class AbstractQuerySet(metaclass=ABCMeta):
    pass
