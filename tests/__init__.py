import unittest
from tests.parsers import ParsingQueries, ParsingExtras
from tests.database import Mongo

import logging

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main(verbosity=1, exit=False)
