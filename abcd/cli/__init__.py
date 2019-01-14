import click


@click.group()
def cli():
    pass


@cli.command('login', short_help='login to the database')
@click.argument('url', default='http://localhost', metavar='<url>')
def login(url):
    """login to the database

    Arguments:\n
    url     The url of abcd api (default: http://localhost)'
    """
    click.echo(f"login: {url}")


@cli.command('list', short_help='list data available in the database')
def db_list():
    """list data available in the database."""
    click.echo(f"list:")


@cli.command('pull', short_help='pull data from the database')
@click.option('--file-format', type=click.Choice(['xyz', 'json']), default='xyz')
def db_pull(file_format):
    """pull data from the database"""
    click.echo(f"push: {file_format}")


@cli.command('push', short_help='push data to the database')
@click.option('-r', '--recursive', is_flag=True)
@click.option('--file-format', type=click.Choice(['xyz', 'json']), default='xyz')
@click.argument('path')
def db_push(recursive, file_format, path):
    """push data to the database

    format
    path      The folder which contains the trajectories
    """
    click.echo(f"push: {recursive}, {file_format}, {path}")


@cli.command('query', short_help='query data from the database')
@click.argument('query_string')
def db_query(query_string):
    """query data from the database"""
    click.echo(f"query: {query_string}")


if __name__ == '__main__':
    cli()

# import argparse
#
# parser = argparse.ArgumentParser(description='Process some integers.')
#
#
# def login_args(args):
#     print(f'login: {args}')
#
#
# def list_args(args):
#     print(f'list: {args}')
#
#
# def pull_args(args):
#     print(f'pull: {args}')
#
#
# def push_args(args):
#     print(f'push: {args}')
#
#
# def query_args(args):
#     print(f'query: {args}')
#
#
# parser.add_argument('-v', '--verbose',
#                     help='Enable verbose mode',
#                     action='store_true')
#
# subparsers = parser.add_subparsers(title='Commands', dest='command')
#
# parser_login = subparsers.add_parser('login', help='login to the database')
# parser_login.set_defaults(func=login_args)
# parser_login.add_argument(dest='url',
#                           help='url of abcd api (default: http://localhost)',
#                           default='http://localhost')
#
# parser_list = subparsers.add_parser('list', help='list data available in the database')
# parser_list.set_defaults(func=list_args)
#
# parser_pull = subparsers.add_parser('pull', help='pull data from the database')
# parser_pull.set_defaults(func=pull_args)
# parser_pull.add_argument(dest='format', choices=['xyz', 'json'])
#
# parser_push = subparsers.add_parser('push', help='pull data from the database')
# parser_push.set_defaults(func=push_args)
# parser_push.add_argument(dest='format', choices=['xyz', 'json'], default='xyz')
# parser_push.add_argument(dest='path', help='folder which contains the trajectories')
# parser_push.add_argument('-r', '--recursive', action='store_true')
#
# parser_query = subparsers.add_parser('query', help='query data from the database')
# parser_query.set_defaults(func=push_args)
# parser_query.add_argument('query_string')
#
#
# def main():
#     args = parser.parse_args()
#
#     if args.command is None:
#         print(parser.format_help())
#         parser.exit()
#
#     if args.verbose:
#         print('verbosity turned on')
#
#     args.func(args)
#
#
# if __name__ == '__main__':
#     main()
