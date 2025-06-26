# Copyright (c) 2025.
# Authors: Tamás K. Stenczel
# This program is distributed under the MIT License, see LICENSE.md.

"""Tests for supporting `mongodb+srv://` URIs"""

from pytest import fixture

from abcd import ABCD


@fixture
def mongo_srv(mocker):
    # mongomock does not pick up what we need, so this is a bespoke mocker
    mock_client = mocker.MagicMock()
    mocker.patch("abcd.backends.atoms_pymongo.MongoClient", mock_client)
    return mock_client


def test_init_mongodb_srv(mongo_srv):
    # client can be created with a mongodb+srv:// URI

    # apparently mongomock breaks if this import is outside
    from abcd.backends.atoms_pymongo import MongoDatabase

    # regression test
    uri = "mongodb+srv://user:pass@democluster.randomstr.mongodb.net/?key=value"

    # create the client
    abcd = ABCD.from_url(uri)

    assert isinstance(abcd, MongoDatabase)
    mongo_srv.assert_called_once_with(host=uri, authSource="admin")
