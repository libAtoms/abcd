import unittest
import logging

from tests.parsers import ParsingQueries, ParsingExtras
from tests.properties import PropertiesTests
from tests.mongo_mock import MongoMock
from tests.opensearch_mock import OpenSearchMock
from tests.opensearch import OpenSearch
from tests.cli import CLI

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    unittest.main(verbosity=1, exit=False)
