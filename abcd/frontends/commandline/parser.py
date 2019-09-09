import logging
from argparse import ArgumentParser, REMAINDER
from abcd.frontends.commandline import commands

logger = logging.getLogger(__name__)

parser = ArgumentParser(description='Command line interface for ABCD database')
parser.add_argument('-v', '--verbose', help='Enable verbose mode', action='store_true')
parser.add_argument('-q', '--query', dest='default_query', action='append', help='Filtering extra quantities',
                    default=[])

parser.add_argument('--read-only', help='Disables all the functions which would modify the database',
                    action='store_true')

subparsers = parser.add_subparsers(title='Commands', dest='command', parser_class=ArgumentParser)

login_parser = subparsers.add_parser('login', help='login to the database')
login_parser.set_defaults(callback_func=commands.login)
login_parser.add_argument('-n', '--name', help='name of the database', default='default')
login_parser.add_argument(dest='url',
                          help='url of abcd api (default: http://localhost)',
                          default='http://localhost')

download_parser = subparsers.add_parser('download', help='download data from the database')
download_parser.set_defaults(callback_func=commands.download)
download_parser.add_argument('-q', '--query', action='append', help='Filtering extra quantities', default=[])
download_parser.add_argument(dest='filename', help='name of the file to store the configurations')

upload_parser = subparsers.add_parser('upload', help='upload any ase supported files to the database')
upload_parser.add_argument('-e', '--extra_infos', action='append', help='Adding extra quantities')
upload_parser.add_argument(dest='path', help='file or folder which contains the xyz files')
upload_parser.add_argument('-i', '--ignore_calc_results', action='store_true',
                           help='Ignore calculators results/parameters')
upload_parser.set_defaults(callback_func=commands.upload)

summary_parser = subparsers.add_parser('summary', help='Discovery mode')
summary_parser.set_defaults(callback_func=commands.summary)
summary_parser.add_argument('-q', '--query', action='append', help='Filtering extra quantities', default=[])
summary_parser.add_argument('-p', '--props', action='append',
                            help='Selecting properties for detailed description')
summary_parser.add_argument('-a', '--all',
                            help='Show everything without truncation of strings and limits of lines',
                            action='store_true', dest='print_all')
summary_parser.add_argument('-n', '--bins', help='The number of bins of the histogram', default=10, type=int)
summary_parser.add_argument('-t', '--trunc',
                            help='Length of string before truncation',
                            default=20, type=int, dest='truncate')

delete_parser = subparsers.add_parser('delete', help='Delete configurations from the database')
delete_parser.set_defaults(callback_func=commands.delete)
delete_parser.add_argument('-q', '--query', action='append', help='Filtering by a query', default=[])
delete_parser.add_argument('-y', '--yes', action='store_true', help='Do the actual deletion.')

key_add_parser = subparsers.add_parser('add-key', help='Adding new key value pairs for a given query')
key_add_parser.set_defaults(callback_func=commands.key_add)
key_add_parser.add_argument('-q', '--query', action='append', help='Filtering by a query', default=[])
key_add_parser.add_argument('-y', '--yes', action='store_true', help='Overwrite?')
key_add_parser.add_argument('keys', help='keys(=value) pairs', nargs='+')

key_rename_parser = subparsers.add_parser('rename-key', help='Rename a specific keys for a given query')
key_rename_parser.set_defaults(callback_func=commands.key_rename)
key_rename_parser.add_argument('-q', '--query', action='append', help='Filtering by a query', default=[])
key_rename_parser.add_argument('-y', '--yes', action='store_true', help='Overwrite?')
key_rename_parser.add_argument('old_keys', help='name of the old key')
key_rename_parser.add_argument('new_keys', help='new name of the key')

key_delete_parser = subparsers.add_parser('delete-key', help='Delete all the keys for a given query')
key_delete_parser.set_defaults(callback_func=commands.key_delete)
key_delete_parser.add_argument('-q', '--query', action='append', help='Filtering by a query', default=[])
key_delete_parser.add_argument('-y', '--yes', action='store_true', help='Do the actual deletion.')
key_delete_parser.add_argument('keys', help='keys(=value) data', nargs='+')

exec_parser = subparsers.add_parser('exec', help='Running custom python code')
exec_parser.set_defaults(callback_func=commands.execute)
exec_parser.add_argument('-q', '--query', action='append', help='Filtering by a query', default=[])
exec_parser.add_argument('-y', '--yes', action='store_true', help='Do the actual execution.')
exec_parser.add_argument('python_code', help='Selecting properties for detailed description')


def main(args=None):
    kwargs = parser.parse_args(args).__dict__

    if kwargs.pop('verbose'):
        # Remove all handlers associated with the root logger object.
        # https://stackoverflow.com/questions/12158048/changing-loggings-basicconfig-which-is-already-set
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO)
        logger.info('Verbose mode is active')

    if not kwargs.pop('command'):
        print(parser.format_help())
        return

    # return namespace.callback_func(namespace)
    callback_func = kwargs.pop('callback_func')
    return callback_func(**kwargs)


if __name__ == '__main__':
    main('summary -q formula~"Si2"'.split())
    main()
    main('upload -e cas -i ../../../tutorials/GB_alphaFe_001/tilt/00110391110_v6bxv2_tv0.4bxv0.2_d1.6z_traj.xyz'.split())
    main('-v login mongodb://mongoadmin:secret@localhost:27017/abcd'.split())
    # main('-v summary'.split())
    # main('-v summary -p energy'.split())
    # main('-v summary -p *'.split())
    # main('add-key -q cas selected user="cas"'.split())
    main('delete-key user'.split())
    # main(['summary', '-p', '*'])
    # main(['summary', '-p', 'info.config_name, info.energy'])
    # main(['summary', '-p', 'info.config_name, info.energy,info.energy;info.energy info.energy'])
    # main(['-s', 'fancy', 'summary', '-p', '*'])
    # main(['summary', '-p', '*'])
    # main(['-s', 'fancy', 'summary'])
    # main(['-v', 'summary', '-p', 'config_type', '-p', 'haha', '-p' 'sdgsdg, dsgsdg,asd fgg', '-q', 'config_type~bcc',
    #      '-q', 'config_type~bcc'])
