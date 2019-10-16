from abcd import ABCD


# from flask_paginate import Pagination, get_page_args


class Database(ABCD):
    """Wrapper for the ABCD factory method for registering a the database for the Flask application."""

    def __init__(self):
        super().__init__()

    def init_app(self, app):
        pass


db = Database()
