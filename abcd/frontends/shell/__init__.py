import logging
import argparse

from abcd import ABCD
from pathlib import Path

from abcd.frontends.shell.config import Config
# from abcd.frontends.shell.parser import ArgumentParser
from abcd.frontends.shell.styles import SimpleStyle, FancyStyle

logger = logging.getLogger(__name__)


def cli():
    return Shell()()


class Shell(argparse.ArgumentParser):

    def __init__(self):

        super().__init__(description='Command line interface for ABCD database')
        self.add_argument('-v', '--verbose', help='Enable verbose mode', action='store_true')
        self.add_argument('-s', '--style', help='style [simple, fancy]', default='simple')

        subparsers = self.add_subparsers(title='Commands', dest='command', parser_class=argparse.ArgumentParser)

        parser_login = subparsers.add_parser('login', help='login to the database')
        parser_login.add_argument('-n', '--name', help='name of the database', default='default')
        parser_login.add_argument(dest='url',
                                  help='url of abcd api (default: http://localhost)',
                                  default='http://localhost')
        parser_login.set_defaults(func=self.login_args)

        parser_download = subparsers.add_parser('download', help='download data from the database')
        parser_download.add_argument(dest='format', choices=['xyz', 'json'])
        parser_download.set_defaults(func=self.download_args)

        parser_upload = subparsers.add_parser('upload', help='upload any ase supported files to the database')
        # parser_upload.add_argument('-r', '--recursive', action='store_true')
        parser_upload.add_argument('-e', '--extra', help='Adding extra quantities', action='store_true')
        parser_upload.add_argument(dest='path', help='file or folder which contains the xyz files')
        parser_upload.set_defaults(func=self.upload_args)

        summary_parser = subparsers.add_parser('summary', help='Discovery mode')
        summary_parser.add_argument('-q', '--query', help='Filtering extra quantities')
        summary_parser.add_argument('-p', '--props', help='Selecting properties for detailed description')
        summary_parser.set_defaults(func=self.summary_args)

        self.config = Config.load()
        # self.parser = ArgumentParser(self.login_args, self.download_args, self.upload_args, self.summary_args)

        # Lazy initialisation of the style and database connection
        self.style = None
        self.db = None

    def init_db(self):
        url = self.config['url']

        if url is None:
            print('Please use abcd login first!')
            self.parser.exit()

        self.db = ABCD(url=url)

    def init_style(self):
        if self.style is None:
            # self.style = SimpleStyle()
            self.style = FancyStyle()

    def login_args(self, args):
        logger.info(f'login args: \n{args}')

        self.db = ABCD(url=args.url)
        self._config.save({
            'url': args.url
        })

    def download_args(self, args):
        self.init_db()

        logger.info(f'download args: \n{args}')

    def upload_args(self, args):
        logger.info(f'upload args: \n{args}')
        self.init_db()

        extra_info = args.extra

        path = Path(args.path)
        if path.is_file():
            self.db.upload()
        if path.is_dir():
            for file in path.glob('.xyz'):
                logger.info(f'Uploaded file: {file}')
                self.db.upload(file, extra_info)
            else:
                logger.info(f'No file found: {path}')
                raise FileNotFoundError()
        else:
            raise FileNotFoundError()

    def summary_args(self, args):
        self.init_db()
        logger.info(f'summary args: \n{args}')

        if args.props is None:
            with self.style as f:

                props = self.db.count_properties(query=args.query)
                f.title('Properties')
                f.h1('Arrays (per atom properties)')
                f.table([[k, v['count']] for k, v in props['arrays'].items()],
                        headers=('label', 'count'),
                        colalign=('left', 'right'),
                        disable_numparse=True,
                        tablefmt='psql')
                f.h1('Infos (properties of the whole configuration)')
                f.table([[k, v['count']] for k, v in props['info'].items()],
                        headers=('label', 'count'),
                        colalign=('left', 'right'),
                        disable_numparse=True,
                        tablefmt='psql')
            return

        elif args.props == '*':
            props = self.db.properties(query=args.query)
            with self.style as f:
                for p in props['arrays']:
                    name = 'arrays.' + p
                    data = self.db.hist(name)
                    f.title(name)
                    if data:
                        f.describe(data)
                        f.hist(data)

                for p in props['info']:
                    name = 'info.' + p
                    data = self.db.hist(name)

                    f.title(name)
                    if data:
                        f.describe(data)
                        f.hist(data)

            return

        else:
            properties = ['info.config_type']

            with self.style as f:
                for p in properties:
                    data = self.db.hist(p)

                    f.title(p)
                    f.hist(data)

            with self.style as f:

                # Count property
                query = {
                    'info.config_type': 'bcc_bulk_54_high'
                }
                query = 'info.config_type=bcc_bulk_54_high'

                f.title('Count query')
                f.print(self.db.count(query))

                # text base histogram
                hist = self.db.hist('info.config_type')
                f.hist(hist)

                data = self.db.hist('info.config_name', query)
                f.hist(data)

                # scalar histogram
                data = self.db.hist('info.energy', query)
                f.hist(data)

    def __call__(self, *args, **kwargs):
        args = self.parse_args(*args, **kwargs)
        if args.verbose:
            logging.basicConfig(level=logging.INFO)
            logger.info('Verbose mode is active')
        else:
            logging.basicConfig(level=logging.ERROR)

        if args.style is 'simple':
            self.style = SimpleStyle()
        elif args.style is 'fancy':
            self.style = FancyStyle()

        if args.command is None:
            print(self.format_help())
        else:
            args.func(args)


if __name__ == '__main__':
    parser = Shell()
    parser()
    parser(['-v', 'summary'])
    # parser(['summary', '-v']) # wrong
    # parser(['summary', '-h']) # ok
    parser(['summary', '-p', '*'])
    # parser(['summary', '-p', 'info.config_name, info.energy'])
    parser(['-s', 'fancy', 'summary', '-p', '*'])
