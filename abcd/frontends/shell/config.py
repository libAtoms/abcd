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
        with open(filename) as json_file:
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
                logger.info(f'Using config file: {file}')
                files.append(file)

            if (Path.cwd() / '.abcd').is_file():
                file = Path.cwd() / '.abcd'
                logger.info(f'Using config file: {file}')
                files.append(file)

            if os.environ.get('ABCD_CONFIG') and Path(os.environ.get('ABCD_CONFIG')).is_file():
                file = Path(os.environ.get('ABCD_CONFIG'))
                logger.info(f'Using config file: {file}')
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

            logger.info(f'Using config file: {file}')

            config = cls.from_json(file)

        return config

    def save(self, local=True):
        file = Path(os.environ.get('ABCD_CONFIG')) if os.environ.get('ABCD_CONFIG') \
            else Path.cwd() / '.abcd' if local \
            else Path.home() / '.abcd'

        logger.info(f'The saved config\'s file: {file}')

        with open(str(file), 'w') as file:
            json.dump(self, file)

    def __repr__(self):
        return f'<{self.__class__.__name__} {dict.__repr__(self)}>'

    def title(self, title, fill='='):
        title = self._trunc(title, self.width - 4)
        template = f'{{:{fill}^20}}'
        return template.format(' ' + title + ' ')

    def table(self, data, hearder=None, style=None):
        template = '{}{}'
        out = '\n'.join(
            template.format(*d) for d in data
        )
        return out

    def histogram(self):
        pass
