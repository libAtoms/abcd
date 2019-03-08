import cmd
import json
import click

from abcd import ABCD
from abcd.config import Config

from ase.io import write


class REPL(cmd.Cmd):
    def __init__(self, ctx):
        cmd.Cmd.__init__(self)
        self.ctx = ctx

    def default(self, line):
        subcommand = cli.commands.get(line)
        if subcommand:
            self.ctx.invoke(subcommand)
        else:
            return cmd.Cmd.default(self, line)


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
    if ctx.invoked_subcommand == 'login':
        return

    config = Config.from_file()

    if config['url'] is None:
        click.echo('Please use abcd login first!')
        raise click.Abort()

    ctx.obj = ABCD(url=config['url'])

    # start a command loop
    if ctx.invoked_subcommand is None:
        repl = REPL(ctx)
        repl.cmdloop()


@cli.command('login', short_help='login to the database')
@click.argument('url', default='http://localhost', metavar='<url>')
@click.argument('collection', default='atoms', metavar='<collection>')
@click.pass_context
def db_login(ctx, url, collection):
    """login to the database

    \b
    Arguments:
    url          The url of abcd api (default: http://localhost)
    collection   The name of collection (default: atoms)
    """
    config = Config.from_file()

    config['url'] = url
    config['collection'] = collection

    # try to connect
    ctx.obj = ABCD(url=config['url'])

    config.save_json()

    # click.echo(f"login: url={url}, collection={collection}")
    click.echo(f"Config file ({config.file.absolute()}) has been updated! ")


@cli.command('list', short_help='list data available in the database')
def db_list():
    """list data available in the database."""
    click.echo(f"list:")


@cli.command('status', short_help='list data available in the database')
@click.pass_obj
def db_status(abcd):
    """list data available in the database."""
    # click.echo(f"connected: {config['url']}")

    abcd.print_info()


@cli.command('properties', short_help='available properties of data in the database')
@click.option('--details', is_flag=True)
@click.option('-q', '--query')
@click.pass_obj
def db_properties(abcd, details, query=None):
    """pull data from the database"""

    if details:
        data = abcd.count_properties(query)
    else:
        data = abcd.properties(query)

    click.echo(json.dumps(data, indent=2))


@cli.command('pull', short_help='pull data from the database')
@click.option('-a', '--append', is_flag=True)
@click.option('-q', '--query')
@click.argument('filename')
@click.pass_obj
def db_pull(abcd, append, query, filename):
    """pull data from the database"""

    atoms = abcd.get_atoms(query)
    write(filename, atoms, append=append)

    click.echo(f"push: filename={filename}")


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


@cli.command('destroy', short_help='destroying the database')
@click.argument('collection')
@click.pass_obj
def db_query(abcd):
    """query data from the database"""
    abcd.destroy()
    click.echo(f"Database destroyed!")


if __name__ == '__main__':
    import sys

    config = Config.from_file()
    config['test1'] = 123
    config['test2'] = 'dsfkjak fds'
    config['test3'] = 'asdfksd kfjs dl'
    config.save_json()

    # sys.argv[1:] = ['status']
    # print(sys.argv)
    # cli()

    # sys.argv[1:] = ['login', 'mongodb://2ef35d3635e9dc5a922a6a42:ac6ce72e259f5ddcc8dd5178@localhost:27017/abcd', 'atoms']
    # print(sys.argv)
    # cli()

    sys.argv[1:] = ['status']
    print(sys.argv)
    cli()

    sys.argv[1:] = ['status']
    print(sys.argv)
    cli()
