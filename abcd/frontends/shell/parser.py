import logging
import argparse

logger = logging.getLogger(__name__)


class ArgumentParser(argparse.ArgumentParser):

    def __init__(self, callback_login, callback_download, callback_upload, callback_summary):
        super().__init__(description='Command line interface for ABCD database')
        self.add_argument('-v', '--verbose', help='Enable verbose mode', action='store_true')
        self.add_argument('-s', '--style', help='style [simple, fancy]', default='simple')
        # self.add_argument('-c', '--show-config', help='Show all the available configs', action='store_true')
        # self.add_argument('-s', '--save_config', help='Generate a local config file', action='store_true')

        subparsers = self.add_subparsers(title='Commands', dest='command', parser_class=argparse.ArgumentParser)

        self.build_login_parser(subparsers, callback_login)
        self.build_download_parser(subparsers, callback_download)
        self.build_upload_parser(subparsers, callback_upload)
        self.build_summary_parser(subparsers, callback_summary)

    @staticmethod
    def build_login_parser(subparsers, callback):
        parser = subparsers.add_parser('login', help='login to the database')
        parser.add_argument('-n', '--name', help='name of the database', default='default')
        parser.add_argument(dest='url',
                            help='url of abcd api (default: http://localhost)',
                            default='http://localhost')
        parser.set_defaults(func=callback)

    @staticmethod
    def build_download_parser(subparsers, callback):
        parser = subparsers.add_parser('download', help='download data from the database')
        parser.add_argument(dest='format', choices=['xyz', 'json'])
        parser.set_defaults(func=callback)

    @staticmethod
    def build_upload_parser(subparsers, callback):
        parser = subparsers.add_parser('upload', help='upload any ase supported files to the database')
        # parser.add_argument('-r', '--recursive', action='store_true')
        parser.add_argument('-e', '--extra', help='Adding extra quantities', action='store_true')
        parser.add_argument(dest='path', help='file or folder which contains the xyz files')
        parser.set_defaults(func=callback)

    @staticmethod
    def build_summary_parser(subparsers, callback):
        parser = subparsers.add_parser('summary', help='Discovery mode')
        parser.add_argument('-q', '--query', help='Filtering extra quantities')
        parser.add_argument('-p', '--props', help='Selecting properties for detailed description')
        parser.set_defaults(func=callback)


if __name__ == '__main__':
    def dummy(args):
        print(args)


    p = ArgumentParser(dummy, dummy, dummy, dummy)
    p.parse_args()
    p.parse_args(['-v', 'login', 'sdfg'])
