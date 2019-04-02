import logging
import argparse
import getpass

from abcd import ABCD
from pathlib import Path

from abcd.frontends.shell.config import Config
from abcd.frontends.shell.styles import SimpleStyle, FancyStyle

logger = logging.getLogger(__name__)


def cli(args=None):
    return ArgumentParser()(args)


class ArgumentParser(argparse.ArgumentParser):

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

        parser_download = subparsers.add_parser('download', help='download data from the database')
        parser_download.add_argument(dest='format', choices=['xyz', 'json'])

        parser_upload = subparsers.add_parser('upload', help='upload any ase supported files to the database')
        # parser_upload.add_argument('-r', '--recursive', action='store_true')
        parser_upload.add_argument('-e', '--extra', help='Adding extra quantities')
        parser_upload.add_argument(dest='path', help='file or folder which contains the xyz files')

        summary_parser = subparsers.add_parser('summary', help='Discovery mode')
        summary_parser.add_argument('-q', '--query', help='Filtering extra quantities')
        summary_parser.add_argument('-p', '--props', help='Selecting properties for detailed description')

    def __call__(self, args=None):

        namespace = self.parse_args(args)

        # Remove all handlers associated with the root logger object.
        # https://stackoverflow.com/questions/12158048/changing-loggings-basicconfig-which-is-already-set
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        if namespace.verbose:
            logging.basicConfig(level=logging.INFO)
            logger.info('Verbose mode is active')
        else:
            logging.basicConfig(level=logging.ERROR)

        if namespace.command is 'login':
            Shell(namespace).login()
        elif namespace.command == 'download':
            Shell(namespace).download()
        elif namespace.command == 'upload':
            Shell(namespace).upload()
        elif namespace.command == 'summary':
            Shell(namespace).summary()
        else:
            print(self.format_help())


class Shell(object):

    def __init__(self, args):
        self.args = args
        self.config = Config.load()

        if args.style == 'simple':
            self.style = SimpleStyle()
        elif args.style == 'fancy':
            self.style = FancyStyle()
        else:
            raise NotImplementedError(f'The style "{args.style}" is not implemented!')

        self.db = None

    def init_db(self):
        url = self.config['url']

        if url is None:
            print('Please use abcd login first!')
            exit()

        self.db = ABCD(url=url)

    def login(self):
        args = self.args
        logger.info(f'login args: \n{args}')

        self.db = ABCD(url=args.url)
        self.config.save({
            'url': args.url
        })

    def download(self):
        args = self.args
        self.init_db()

        logger.info(f'download args: \n{args}')

    def upload(self):
        args = self.args
        logger.info(f'upload args: \n{args}')
        self.init_db()

        extra_info = args.extra

        path = Path(args.path)
        if path.is_file():
            self.db.upload(path, extra_info)
        elif path.is_dir():
            for file in path.glob('.xyz'):
                logger.info(f'Uploaded file: {file}')
                self.db.upload(file, extra_info)
            else:
                logger.info(f'No file found: {path}')
                raise FileNotFoundError()
        else:
            raise FileNotFoundError()

    def summary(self):
        self.init_db()

        args = self.args
        logger.info(f'summary args: \n{args}')

        if args.props is None:
            with self.style as f:

                props = self.db.count_properties(query=args.query)
                if props['arrays']:
                    f.title('Properties')
                    f.h1('Arrays (per atom properties)')

                    labels, counts = [], []
                    for k, v in props['arrays'].items():
                        labels.append(k)
                        counts.append(v['count'])

                    f.hist({
                        'type': 'hist_labels',
                        'labels': labels,
                        'counts': counts
                    })

                if props['info']:
                    f.h1('Infos (properties of the whole configuration)')

                    labels, counts = [], []
                    for k, v in props['info'].items():
                        labels.append(k)
                        counts.append(v['count'])

                    f.hist({
                        'type': 'hist_labels',
                        'labels': labels,
                        'counts': counts
                    })
            return

        elif args.props == '*':
            props = self.db.properties(query=args.query)
            with self.style as f:
                for p in props['arrays']:
                    name = 'arrays.' + p
                    data = self.db.hist(name, query=args.query)
                    f.title(name)
                    if data:
                        f.describe(data)
                        f.hist(data)

                for p in props['info']:
                    name = 'info.' + p
                    data = self.db.hist(name, query=args.query)

                    f.title(name)
                    if data:
                        f.describe(data)
                        f.hist(data)

            return

        def parse_properties(s):
            import re
            return re.split(r';\s*|,\s*|\s+', s)

        properties = parse_properties(args.props)

        with self.style as f:
            for p in properties:

                data = self.db.hist('arrays.' + p, query=args.query)

                if data:
                    f.title('arrays.' + p)
                    f.describe(data)
                    f.hist(data)

                data = self.db.hist('info.' + p, query=args.query)

                if data:
                    f.title('info.' + p)
                    f.describe(data)
                    f.hist(data)


if __name__ == '__main__':
    cli()
    cli(['summary'])
    cli(['-v', 'summary'])
    # cli(['summary', '-v'])  # wrong
    # cli(['summary', '-h'])  # ok
    # cli(['summary', '-p', '*'])
    cli(['summary', '-p', 'info.config_name, info.energy'])
    cli(['summary', '-p', 'info.config_name, info.energy,info.energy;info.energy info.energy'])
    cli(['-s', 'fancy', 'summary', '-p', '*'])
    cli(['summary', '-p', '*'])
    cli(['-s', 'fancy', 'summary'])
