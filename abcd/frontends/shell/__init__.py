import logging

from abcd import ABCD
from pathlib import Path

from abcd.frontends.shell.config import Config
from abcd.frontends.shell.parser import ArgumentParser
from abcd.frontends.shell.styles import SimpleStyle

logger = logging.getLogger(__name__)


def cli():
    return Shell()()


class Shell(object):

    def __init__(self):

        self.config = Config.load()
        self.format = SimpleStyle()
        self.parser = ArgumentParser(self.login_args, self.download_args, self.upload_args, self.summary_args)

        self._db = None  # Lazy initialisation of the config and database connection

    def abcd_init(self):
        url = self.config['url']

        if url is None:
            print('Please use abcd login first!')
            self.parser.exit()

        self._db = ABCD(url=url)

    def login_args(self, args):
        logger.info(f'login args: \n{args}')

        self._db = ABCD(url=args.url)
        self._config.save({
            'url': args.url
        })

    def download_args(self, args):
        self.abcd_init()

        logger.info(f'download args: \n{args}')

    def upload_args(self, args):
        logger.info(f'upload args: \n{args}')
        self.abcd_init()

        extra_info = args.extra

        path = Path(args.path)
        if path.is_file():
            self._db.upload()
        if path.is_dir():
            for file in path.glob('.xyz'):
                logger.info(f'Uploaded file: {file}')
                self._db.upload(file, extra_info)
            else:
                logger.info(f'No file found: {path}')
                raise FileNotFoundError()
        else:
            raise FileNotFoundError()

    def summary_args(self, args):
        self.abcd_init()
        logger.info(f'summary args: \n{args}')

        if args.props is None:
            props = self._db.properties()
            properties = list(f'info.{p}' for p in props['info'])

        else:
            properties = ['info.config_type']

        with self.format as f:
            for p in properties:
                data = self._db.property(p)

                f.title(p)
                f.hist(data)

        with self.format as f:
            props = self._db.properties()
            f.title('Properties')
            f.h1('Arrays (per atom properties)')
            f.print(props['arrays'])
            f.h1('Infos (properties of the hole configuration)')
            f.print(props['info'])

            props = self._db.count_properties()
            f.title('Properties')
            f.h1('Arrays (per atom properties)')
            f.table([[k, v['count']] for k, v in props['arrays'].items()],
                    headers=('label', 'count'),
                    colalign=('left', 'right'),
                    disable_numparse=True,
                    tablefmt='psql')
            f.h1('Infos (properties of the hole configuration)')
            f.table([[k, v['count']] for k, v in props['info'].items()],
                    headers=('label', 'count'),
                    colalign=('left', 'right'),
                    disable_numparse=True,
                    tablefmt='psql')

            # Count property
            query = {
                'info.config_type': 'bcc_bulk_54_high'
            }
            query = 'info.config_type=bcc_bulk_54_high'

            f.title('Count query')
            f.print(self._db.count(query))

            # text base histogram
            data = self._db.property('info.config_type')
            f.title('Histogram: info.config_type')
            f.h1('Property: info.config_type')
            f.hist(data)

            data = self._db.property('info.config_name', query)
            f.title('Histogram: info.config_name')
            f.h1('Property: info.config_name')
            f.hist(data)

            # scalar histogram
            data = self._db.property('info.energy', query)
            f.title('Histogram: info.energy')
            f.h1('Property: info.energy')
            f.hist(data)

    def __call__(self, *args, **kwargs):
        self.parser.parse_args(*args, **kwargs)


if __name__ == '__main__':
    parser = Shell()
    parser(['-v', 'summary'])

    cli()
