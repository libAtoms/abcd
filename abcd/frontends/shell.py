import os
import json
import shutil
import logging
import argparse
import numpy as np

from abcd import ABCD
from pathlib import Path
from pprint import pprint
from collections import Counter
from abc import ABCMeta, abstractmethod

from tabulate import tabulate

logger = logging.getLogger(__name__)


class Parser(object):
    def properties(self, s):
        # * or comma/space/; separated list of identifiers
        pass

    def extra(self, s):
        # key=value pairs for extra data for info
        # simplified version of querystring
        pass

    def querystring(self, s):
        pass


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
                raise FileNotFoundError()

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


class Format(object):
    partialBlocks = ["▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"]  # char=pb

    def __init__(self, width=None, height=None):
        self.width = width
        self.height = height

    def __enter__(self):
        if self.width is None:
            self.width = shutil.get_terminal_size()[0]

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def title(self, title):
        fill = '='
        template = f'{{:{fill}^{self.width}}}'
        title = self._trunc(title, self.width - 4)
        print('', template.format(' ' + title + ' '), sep=os.linesep)

    def h1(self, title):
        title = self._trunc(title, self.width)
        print('', title, '=' * len(title), sep=os.linesep)

    def h2(self, title):
        title = self._trunc(title, self.width)
        print('', title, '-' * len(title), sep=os.linesep)

    def line(self, data):
        print(data)

    def list(self, data):
        print(*(' - ' + item for item in data), sep=os.linesep)

    def hist(self, data, **kwargs):
        if isinstance(data, list):
            if isinstance(data[0], float):
                self.hist_float(data, **kwargs)
            elif isinstance(data[0], str):
                self.hist_str(data, **kwargs)
            else:
                raise NotImplementedError(f'Histogram for list of {type(data)} types are not supported!')
        else:
            raise NotImplementedError(f'Histogram for {type(data)} types are not supported!')

    def hist_str(self, data, max_width=80):
        data = Counter(data)

        keys = data.keys()
        values = data.values()

        self.h2('Summary')
        self.table(
            ((f'{sum(values)}', f'{len(keys)}', f'{min(values)}', f'{max(values)}', f'{sum(values) / len(keys):.8g}'),),
            headers=('total', '# of categories', 'min', 'max', 'mean'), disable_numparse=True, tablefmt="grid")

        self.h2('Histogram')

        def get_width_hack(keys, values):
            # generating a table just to measure its width
            return len(tabulate(
                list([item, f'{count}'] for item, count in zip(keys, values)),
                headers=('Name', 'Count'),
                colalign=("left", "right"),
                disable_numparse=True,
                tablefmt="psql"
            ).split('\n', 1)[0])

        table_width = get_width_hack(keys, values)

        if table_width > self.width:
            self.print('Too small terminal!')
        else:
            max_width = self.width - table_width - 3

            ratio = max_width / max(values)
            scales = (int(ratio * value) for value in values)

            self.table(
                list([item, f'{count}', '▉' * scale] for item, count, scale in zip(keys, values, scales)),
                headers=('Name', 'Count', 'Histogram'),
                colalign=("left", "right", "left"),
                disable_numparse=True,
                tablefmt="psql"
            )

    def hist_float(self, data, bins=10):
        data = np.array(data)

        self.h2('Summary')
        self.table(
            ((f'{len(data)}', f'{data.min():.8g}', f'{data.max():.8g}', f'{data.mean():.8g}', f'{data.std():.8g}',
              f'{data.var():.8g}'),),
            headers=('total', 'min', 'max', 'mean', 'std', 'var'), disable_numparse=True, tablefmt="grid")

        self.h2('Histogram')

        def get_width_hack(hist, bin_edges):
            # generating a table just to measure its width
            return len(tabulate(
                list([f'{lower:.8g}', f'{upper:.8g}', f'{count}']
                     for lower, upper, count in zip(bin_edges[:-1], bin_edges[1:], hist)),
                headers=('Lower', 'Upper', 'Count'),
                colalign=("decimal", "decimal", "right"),
                disable_numparse=True,
                tablefmt="psql"
            ).split('\n', 1)[0])

        hist, bin_edges = np.histogram(data, bins=bins)
        table_width = get_width_hack(hist, bin_edges)

        if table_width > self.width:
            self.print('Too small terminal!')
        else:
            max_width = self.width - table_width - 3
            scales = (max_width / hist.max() * hist).astype(int)

            self.table(
                list([f'{lower:.8g}', f'{upper:.8g}', f'{count}', '▉' * scale]
                     for lower, upper, count, scale in zip(bin_edges[:-1], bin_edges[1:], hist, scales)),
                headers=('Lower', 'Upper', 'Count', 'Histogram'),
                colalign=("decimal", "decimal", "right", "left"),
                disable_numparse=True,
                tablefmt="psql"
            )

    @staticmethod
    def print(*args, **kwargs):
        print(*args, **kwargs)

    def pprint(self, *args):
        pprint(*args, width=self.width)

    @staticmethod
    def table(*args, **kwargs):
        print(tabulate(*args, **kwargs))

    def _trunc(self, text, width=None):
        width = width if width else self.width
        return text if len(text) < width else text[:width - 3] + '...'

    @staticmethod
    def newline():
        print('')


class ArgumentParser(metaclass=ABCMeta):

    def __init__(self):
        self.parser = argparse.ArgumentParser(description='Command line interface for ABCD database')
        self.parser.add_argument('-v', '--verbose', help='Enable verbose mode', action='store_true')
        self.parser.add_argument('-c', '--show-config', help='Show all the available configs', action='store_true')
        self.parser.add_argument('-s', '--save_config', help='Generate a local config file', action='store_true')

        self.subparsers = self.parser.add_subparsers(title='Commands', dest='command')

        self.build_login_parser()
        self.build_download_parser()
        self.build_upload_parser()
        self.build_summary_parser()

        self._db = None
        self._config = None

    def build_login_parser(self):

        parser = self.subparsers.add_parser('login', help='login to the database')
        parser.add_argument('-n', '--name', help='name of the database', default='default')
        parser.add_argument(dest='url',
                            help='url of abcd api (default: http://localhost)',
                            default='http://localhost')
        parser.set_defaults(func=self.login_args)

    def build_download_parser(self):
        parser = self.subparsers.add_parser('download', help='download data from the database')
        parser.add_argument(dest='format', choices=['xyz', 'json'])
        parser.set_defaults(func=self.download_args)

    def build_upload_parser(self):
        parser = self.subparsers.add_parser('upload', help='upload any ase supported files to the database')
        # parser.add_argument('-r', '--recursive', action='store_true')
        parser.add_argument('-e', '--extra', help='Adding extra quantities', action='store_true')
        parser.add_argument(dest='path', help='file or folder which contains the xyz files')
        parser.set_defaults(func=self.upload_args)

    def build_summary_parser(self):
        parser = self.subparsers.add_parser('summary', help='Discovery mode')
        parser.add_argument('-q', '--query', help='Filtering extra quantities')
        parser.add_argument('-p', '--props', help='Selecting properties for detailed description')
        parser.set_defaults(func=self.summary_args)

    def main(self):
        args = self.parser.parse_args()

        if args.command is None:
            print(self.parser.format_help())
            if args.verbose:
                self.all_help()
            self.parser.exit()

        if args.verbose:
            logging.basicConfig(level=logging.INFO)
            logger.info('Verbose mode is active')

        args.func(args)

    def all_help(self):
        """only for debugging"""
        parser = self.parser

        # retrieve subparsers from parser
        subparsers_actions = [
            action for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)]

        # there will probably only be one subparser_action,
        # but better save than sorry
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                print(f"{f' Sub-parser {choice} ':=^60}")
                print(subparser.format_help())

    @abstractmethod
    def login_args(self, args):
        pass

    @abstractmethod
    def download_args(self, args):
        pass

    @abstractmethod
    def upload_args(self, args):
        pass

    @abstractmethod
    def summary_args(self, args):
        pass


class Shell(ArgumentParser):

    def __init__(self):
        super().__init__()

        # Lazy initialisation of the config and database connection
        self._db = None
        self._config = None

    def abcd_init(self):
        config = Config.load()
        url = config['url']

        if url is None:
            print('Please use abcd login first!')
            self.parser.exit()
            # exit()

        self._config = config
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

        with Format() as f:
            props = self._db.properties()
            f.title('Properties')
            f.h1('Arrays (per atom properties)')
            f.list(props['arrays'])
            f.h1('Infos (properties of the hole configuration)')
            f.list(props['info'])

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


def cli():
    parser = Shell()
    return parser.main()


if __name__ == '__main__':
    parser = Shell()
    parser.main()

# TODO: REPL
# class Cmd(cmd.Cmd):
#
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.error = False
#
#     def do_select(self, arg):
#         print(self.select(arg))
#
#     def do_count(self, arg):
#         print(self.count(arg))
#
#     def default(self, arg):
#         super().default(arg)
#         self.error = True
#
#
# def repl():
#     repl = Cmd()
#     repl.cmdloop()
