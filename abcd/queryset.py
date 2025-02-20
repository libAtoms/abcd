from abc import ABCMeta
import logging

logger = logging.getLogger(__name__)


class Query(dict):
    pass


class AbstractQuerySet(metaclass=ABCMeta):  # noqa: B024
    pass
