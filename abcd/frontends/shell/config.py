import os
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class Config(dict):
    # https://martin-thoma.com/configuration-files-in-python/
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_json(cls, filename):
        with open(str(filename)) as json_file:
            obj = json.loads(json_file.read())

        return cls(obj)

    @classmethod
    def load(cls, cumulative=True):

        # TODO: Precedences: Should it be cumulative(updated) of exclusive?
        if cumulative:
            config = cls()
            files = []

            # Load all the defaults
            if (Path.home() / '.abcd').is_file():
                file = Path.home() / '.abcd'
                logger.info('Using config file: {}'.format(file))
                files.append(file)

            if (Path.cwd() / '.abcd').is_file():
                file = Path.cwd() / '.abcd'
                logger.info('Using config file: {}'.format(file))
                files.append(file)

            if os.environ.get('ABCD_CONFIG') and Path(os.environ.get('ABCD_CONFIG')).is_file():
                file = Path(os.environ.get('ABCD_CONFIG'))
                logger.info('Using config file: {}'.format(file))
                files.append(file)

            if files:
                for file in files:
                    config.update(cls.from_json(file))
            else:
                logger.warning('No config found!')
                # raise FileNotFoundError()

        else:
            if os.environ.get('ABCD_CONFIG') and Path(os.environ.get('ABCD_CONFIG')).is_file():
                file = Path(os.environ.get('ABCD_CONFIG'))
            elif (Path.cwd() / '.abcd').is_file():
                file = Path.cwd() / '.abcd'
            elif (Path.home() / '.abcd').is_file():
                file = Path.home() / '.abcd'
            else:
                raise FileNotFoundError()

            logger.info('Using config file: {}'.format(file))

            config = cls.from_json(file)

        return config

    def save(self, local=True):
        file = Path(os.environ.get('ABCD_CONFIG')) if os.environ.get('ABCD_CONFIG') \
            else Path.cwd() / '.abcd' if local \
            else Path.home() / '.abcd'

        logger.info('The saved config\'s file: {}'.format(file))

        with open(str(file), 'w') as file:
            json.dump(self, file)

    def __repr__(self):
        return '<{} {}>'.format(self.__class__.__name__, dict.__repr__(self))
