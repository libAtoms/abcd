import logging

from abcd import ABCD

from abcd.frontends.cli.config import Config
from abcd.frontends.cli.styles import SimpleStyle, FancyStyle
from abcd.parsers.queries import QueryParser

logger = logging.getLogger(__name__)
parser = QueryParser()


def init_style(func):
    def wrapper(*args, **kwargs):

        style_str = kwargs.pop('style')
        if style_str == 'simple':
            style = SimpleStyle()
        elif style_str == 'fancy':
            style = FancyStyle()
        else:
            print('The style "{}" is not implemented!'.format(style_str))
            exit(1)

        func(*args, style=style, **kwargs)

    return wrapper


def init_config(func):
    config = Config.load()

    def wrapper(*args, **kwargs):
        func(*args, config=config, **kwargs)

    return wrapper


def init_db(func):
    def wrapper(*args, config, **kwargs):
        url = config.get('url', None)

        if url is None:
            print('Please use abcd login first!')
            exit(1)

        db = ABCD.from_url(url=url)

        # TODO: clean this part
        # TODO: resturn with a new object (ASTDict)
        default_query = kwargs.pop('default_query', [])
        query = kwargs.pop('query', [])
        q = default_query if default_query else []
        q += query if query else []

        q = parser.parse(q)

        func(*args, db=db, query=q, **kwargs)

    return wrapper


def check_readonly(func):
    def wrapper(*args, **kwargs):
        read_only = kwargs.pop('read_only')
        if read_only:
            print('In read only mode, you can\'t modify the data in the database')
            exit(1)

        func(*args, **kwargs)

    return wrapper
