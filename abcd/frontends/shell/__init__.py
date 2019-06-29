import re
import logging
import argparse

from abcd import ABCD
from pathlib import Path

from abcd.frontends.shell.config import Config
from abcd.frontends.shell.styles import SimpleStyle, FancyStyle

from abcd.backends.abstract import URLError, AuthenticationError

logger = logging.getLogger(__name__)


def cli(args=None):
    return ArgumentParser()(args)


class ArgumentParser(argparse.ArgumentParser):

    def __init__(self):
        super().__init__(description='Command line interface for ABCD database')
        self.add_argument('-v', '--verbose', help='Enable verbose mode', action='store_true')
        # self.add_argument('--remote-access-query')  no upload, no cache
        self.add_argument('-s', '--style', help='style [simple, fancy]', default='simple')
        self.add_argument('-q', '--query', dest='default_query', action='append', help='Filtering extra quantities')

        subparsers = self.add_subparsers(title='Commands', dest='command', parser_class=argparse.ArgumentParser)

        parser_login = subparsers.add_parser('login', help='login to the database')
        parser_login.add_argument('-n', '--name', help='name of the database', default='default')
        parser_login.add_argument(dest='url',
                                  help='url of abcd api (default: http://localhost)',
                                  default='http://localhost')

        parser_download = subparsers.add_parser('download', help='download data from the database')
        parser_download.add_argument('-q', '--query', action='append', help='Filtering extra quantities')
        # parser_download.add_argument(dest='format', choices=['xyz', 'json'], default='xyz')
        parser_download.add_argument(dest='filename', help='name of the file to store the configurations')

        parser_upload = subparsers.add_parser('upload', help='upload any ase supported files to the database')
        # parser_upload.add_argument('-r', '--recursive', action='store_true')
        parser_upload.add_argument('-e', '--extra', action='append', help='Adding extra quantities')
        parser_upload.add_argument(dest='path', help='file or folder which contains the xyz files')

        summary_parser = subparsers.add_parser('summary', help='Discovery mode')
        summary_parser.add_argument('-q', '--query', action='append', help='Filtering extra quantities')
        summary_parser.add_argument('-p', '--props', action='append',
                                    help='Selecting properties for detailed description')
        summary_parser.add_argument('-a', '--all',
                                    help='Show everything without truncation of strings and limits of lines',
                                    action='store_true')
        summary_parser.add_argument('-n', '--bins', help='The number of bins of the histogram', default=10, type=int)
        summary_parser.add_argument('-t', '--trunc', help='Length of string before truncation', default=20, type=int)

        delete_parser = subparsers.add_parser('delete', help='Delete configurations from the database')
        delete_parser.add_argument('-q', '--query', action='append', help='Filtering by a query')
        delete_parser.add_argument('-y', '--yes', action='store_true', help='Do the actual deletion.')

        key_add_parser = subparsers.add_parser('key-add', help='Adding new key value pairs for a given query')
        key_add_parser.add_argument('-q', '--query', action='append', help='Filtering by a query')
        key_add_parser.add_argument('-y', '--yes', action='store_true', help='Overwrite?')
        key_add_parser.add_argument('keys', help='keys(=value) pairs', nargs='+')

        key_rename_parser = subparsers.add_parser('key-rename', help='Rename a specific keys for a given query')
        key_rename_parser.add_argument('-q', '--query', action='append', help='Filtering by a query')
        key_rename_parser.add_argument('-y', '--yes', action='store_true', help='Overwrite?')
        key_rename_parser.add_argument('old_keys', help='name of the old key')
        key_rename_parser.add_argument('new_keys', help='new name of the key')

        key_delete_parser = subparsers.add_parser('key-delete', help='Delete all the keys for a given query')
        key_delete_parser.add_argument('-q', '--query', action='append', help='Filtering by a query')
        key_delete_parser.add_argument('-y', '--yes', action='store_true', help='Do the actual deletion.')
        key_delete_parser.add_argument('keys', help='keys(=value) data', nargs='+')

        exec_parser = subparsers.add_parser('exec', help='Running custom python code')
        exec_parser.add_argument('-q', '--query', action='append', help='Filtering by a query')
        exec_parser.add_argument('-y', '--yes', action='store_true', help='Do the actual execution.')
        exec_parser.add_argument('python_code', help='Selecting properties for detailed description')

        cache_parser = subparsers.add_parser('cache', help='Caching data from remote databases')
        cache_parser.add_argument('-q', '--query', action='append', help='Filtering by a query')
        cache_parser.add_argument('-p', '--props', help='Selecting properties for detailed description')

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

        if namespace.command == 'login':
            Shell(namespace).login()
        elif namespace.command == 'download':
            Shell(namespace).download()
        elif namespace.command == 'upload':
            Shell(namespace).upload()
        elif namespace.command == 'summary':
            Shell(namespace).summary()
        elif namespace.command == 'tag':
            Shell(namespace).tag()
        elif namespace.command == 'delete':
            Shell(namespace).delete()
        elif namespace.command == 'key-add':
            Shell(namespace).key_add()
        elif namespace.command == 'key-rename':
            Shell(namespace).key_rename()
        elif namespace.command == 'key-delete':
            Shell(namespace).key_delete()
        elif namespace.command == 'exec':
            Shell(namespace).exec()
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
            raise NotImplementedError('The style "{}" is not implemented!'.format(args.style))

        if args.command == 'login':
            self.db = None
        else:
            url = self.config.get('url', None)

            if url is None:
                print('Please use abcd login first!')
                exit()
            self.db = ABCD(url=url)

    def login(self):
        args = self.args
        logger.info('login args: \n{}'.format(args))

        try:
            self.db = ABCD(url=args.url)
        except URLError:
            print('Wrong connection: Please check the parameters of th url!')
            exit(1)
        except AuthenticationError:
            print('Authentication failed.')
            exit(1)

        print('Successfully connected to the database!')
        print(" type:       {type}\n"
              " hostname:   {host}\n"
              " port:       {port}\n"
              " database:   {db}\n"
              " # of confs: {number of confs}".format(**self.db.info()))

        self.config['url'] = args.url
        self.config.save()

    def download(self):
        args = self.args
        logger.info('download args: \n{}'.format(args))

        from ase.io import write
        filename = args.filename

        query = []
        query += args.default_query if args.default_query else []
        query += args.query if args.query else []

        write(filename, list(self.db.get_atoms(query=query)))

    def delete(self):
        args = self.args
        logger.info('delete args: \n{}'.format(args))

        query = []
        query += args.default_query if args.default_query else []
        query += args.query if args.query else []

        if not args.yes:
            print('Please use --yes for deleting {} configurations'.format(self.db.count(query=query)))
            exit()

        count = self.db.delete(query=query)
        print('{} configuration has been deleted'.format(count))

    def upload(self):
        args = self.args
        logger.info('upload args: \n{}'.format(args))

        extra_info = args.extra

        path = Path(args.path)
        if path.is_file():
            self.db.upload(path, extra_info)
        elif path.is_dir():
            for file in path.glob('.xyz'):
                logger.info('Uploaded file: {}'.format(file))
                self.db.upload(file, extra_info)
            else:
                logger.info('No file found: {}'.format(path))
                raise FileNotFoundError()
        else:
            raise FileNotFoundError()

    def summary(self):
        args = self.args

        query = []
        query += args.default_query if args.default_query else []
        query += args.query if args.query else []

        logger.info('summary args: \n{}'.format(self.args))

        if self.args.all:
            bins, truncate = None, None
        else:
            bins, truncate = self.args.bins, self.args.trunc

        with self.style as f:

            if self.args.props is None:
                props_list = None
            else:
                props_list = []
                for prop in self.args.props:
                    props_list.extend(re.split(r';\s*|,\s*|\s+', prop))

                if '*' in props_list:
                    props_list = '*'

                logging.info('property list: {}'.format(props_list))

            if props_list is None:

                total = self.db.count(query)
                print('Total number of configurations: {}'.format(total))

                props = self.db.count_properties(query=query)

                if props['arrays']:
                    f.title('Properties')
                    f.h1('Arrays (per atom properties)')

                    labels, counts = [], []
                    for k in sorted(props['arrays'], key=str.lower):
                        labels.append(k)
                        counts.append(props['arrays'][k]['count'])

                    f.hist({
                        'type': 'hist_labels',
                        'labels': labels,
                        'counts': counts
                    })

                if props['info']:
                    f.h1('Infos (properties of the whole configuration)')

                    labels, counts = [], []
                    for k in sorted(props['info'], key=str.lower):
                        labels.append(k)
                        counts.append(props['info'][k]['count'])

                    f.hist({
                        'type': 'hist_labels',
                        'labels': labels,
                        'counts': counts
                    })

                if props['derived']:
                    f.h1('Derived')

                    labels, counts = [], []
                    for k in sorted(props['derived'], key=str.lower):
                        labels.append(k)
                        counts.append(props['derived'][k]['count'])

                    f.hist({
                        'type': 'hist_labels',
                        'labels': labels,
                        'counts': counts
                    })

            elif props_list == '*':
                props = self.db.properties(query=query)

                for p in props['arrays']:
                    name = 'arrays.' + p
                    data = self.db.hist(name, query=query, bins=bins, truncate=truncate)
                    f.title(name)
                    if data:
                        f.describe(data)
                        f.hist(data)

                for p in props['info']:
                    name = 'info.' + p
                    data = self.db.hist(name, query=query, bins=bins, truncate=truncate)

                    f.title(name)
                    if data:
                        f.describe(data)
                        f.hist(data)
            else:
                for p in props_list:
                    data = self.db.hist('arrays.' + p, query=query, bins=bins, truncate=truncate)

                    if data:
                        f.title('arrays.' + p)
                        f.describe(data)
                        f.hist(data)

                    data = self.db.hist('info.' + p, query=query, bins=bins, truncate=truncate)

                    if data:
                        f.title('info.' + p)
                        f.describe(data)
                        f.hist(data)

                    data = self.db.hist('derived.' + p, query=query, bins=bins, truncate=truncate)

                    if data:
                        f.title('derived.' + p)
                        f.describe(data)
                        f.hist(data)

    def key_add(self):
        args = self.args
        logger.info('tag args: \n{}'.format(self.args))

    def key_rename(self):
        args = self.args
        logger.info('tag args: \n{}'.format(self.args))

    def key_delete(self):
        args = self.args
        logger.info('tag args: \n{}'.format(self.args))

    def exec(self):
        args = self.args
        logger.info('tag args: \n{}'.format(self.args))

        query = []
        query += args.default_query if args.default_query else []
        query += args.query if args.query else []

        if not args.yes:
            print('Please use --yes for executing code on {} configurations'.format(self.db.count(query=query)))
            exit()

        self.db.exec(args.python_code, query)


if __name__ == '__main__':
    cli()
    # cli(['-v', 'summary', '-q', 'username=fekad virial dummy'])
    # cli(['-v', 'summary', '-q', 'username=fekad virial dummy', '-q', 'virial dummy'])
    # exit()
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
    cli(['-v', 'summary', '-p', 'config_type', '-p', 'haha', '-p' 'sdgsdg, dsgsdg,asd fgg', '-q', 'config_type~bcc',
         '-q', 'config_type~bcc'])
