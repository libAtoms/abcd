import logging

import os
import numpy as np

from abcd.frontends.commandline.decorators import init_db, init_config, check_readonly

logger = logging.getLogger(__name__)


@init_config
def login(*, config, name, url, **kwargs):
    logger.info('login args: \nconfig:{}, name:{}, url:{}, kwargs:{}'.format(config, name, url, kwargs))
    from abcd import ABCD
    from abcd.errors import URLError, AuthenticationError

    try:
        db = ABCD.from_url(url=url)
        info = db.info()

        config['url'] = url
        config.save()

        print('Successfully connected to the database!')
        print(" type:       {type}\n"
              " hostname:   {host}\n"
              " port:       {port}\n"
              " database:   {db}\n"
              " # of confs: {number of confs}".format(**info))

    except URLError:
        print('Wrong connection: Please check the parameters of th url!')
        exit(1)
    except AuthenticationError:
        print('Authentication failed!')
        exit(1)


@init_config
@init_db
def download(*, db, query, filename, **kwargs):
    logger.info('download\n kwargs: {}'.format(kwargs))

    from ase.io import write
    write(filename, list(db.get_atoms(query=query)))


@init_config
@init_db
@check_readonly
def delete(*, db, query, yes, **kwargs):
    logger.info('delete\n kwargs: {}'.format(kwargs))

    if not yes:
        print('Please use --yes for deleting {} configurations'.format(db.count(query=query)))
        exit(1)

    count = db.delete(query=query)
    print('{} configuration has been deleted'.format(count))


@init_config
@init_db
@check_readonly
def upload(*, db, path, extra_info, **kwargs):
    from pathlib import Path

    path = Path(path)
    if path.is_file():
        db.upload(path, extra_info)
    elif path.is_dir():
        for file in path.glob('.xyz'):
            logger.info('Uploaded file: {}'.format(file))
            db.upload(file, extra_info)
        else:
            logger.info('No file found: {}'.format(path))
            raise FileNotFoundError()
    else:
        raise FileNotFoundError()


@init_config
@init_db
def summary(*, db, query, print_all, bins, truncate, props, **kwargs):
    logger.info('summary\n kwargs: {}'.format(kwargs))
    logger.info('query: {}'.format(query))

    if print_all:
        bins, truncate = None, None

    if props is None:
        props_list = None
    else:
        import re

        props_list = []
        for prop in props:
            # TODO: Check that is this the right place?
            props_list.extend(re.split(r';\s*|,\s*|\s+', prop))

        if '*' in props_list:
            props_list = '*'

        logging.info('property list: {}'.format(props_list))

    f = SimpleStyle()
    if props_list is None:

        total = db.count(query)
        props = db.count_properties(query=query)

        print('Total number of configurations: {}'.format(total))

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
        props = db.properties(query=query)

        for p in props['arrays']:
            name = 'arrays.' + p
            data = db.hist(name, query=query, bins=bins, truncate=truncate)

            f.title(name)
            if data:
                f.describe(data)
                f.hist(data)

        for p in props['info']:
            name = 'info.' + p
            data = db.hist(name, query=query, bins=bins, truncate=truncate)

            f.title(name)
            if data:
                f.describe(data)
                f.hist(data)

    else:
        for p in props_list:
            data = db.hist('arrays.' + p, query=query, bins=bins, truncate=truncate)

            if data:
                f.title('arrays.' + p)
                f.describe(data)
                f.hist(data)

            data = db.hist('info.' + p, query=query, bins=bins, truncate=truncate)

            if data:
                f.title('info.' + p)
                f.describe(data)
                f.hist(data)

            data = db.hist('derived.' + p, query=query, bins=bins, truncate=truncate)

            if data:
                f.title('derived.' + p)
                f.describe(data)
                f.hist(data)


@check_readonly
@init_config
@init_db
def key_add(*, db, query, keys, **kwargs):
    from abcd.parsers.extras_fallback import key_val_str_to_dict

    data = key_val_str_to_dict(' '.join(keys))

    if db.count(query=query + [' or '.join(data.keys())]):
        print('The new key already exist for the given query! '
              'Please make sure that the target key name don\'t exist')
        exit(1)

    db.add_property(data, query=query)


@check_readonly
@init_config
@init_db
def key_rename(*, db, query, old_keys, new_keys, **kwargs):
    if db.count(query=query + [old_keys, new_keys]):
        print('The new key already exist for the given query! '
              'Please make sure that the target key name don\'t exist')
        exit(1)

    db.rename_property(old_keys, new_keys, query=query)


@check_readonly
@init_config
@init_db
def key_delete(*, db, query, yes, keys, **kwargs):
    query += keys

    if not yes:
        print('Please use --yes for deleting keys from {} configurations'.format(db.count(query=query)))
        exit(1)

    for k in keys:
        db.delete_property(k, query=query)


@check_readonly
@init_config
@init_db
def execute(*, db, query, yes, python_code, **kwargs):
    if not yes:
        print('Please use --yes for executing code on {} configurations'.format(db.count(query=query)))
        exit(1)

    db.exec(python_code, query)


class SimpleStyle(object):
    partialBlocks = ["▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"]  # char=pb

    def __init__(self, width=80):
        self.width = width

    def title(self, title):
        # template = '{{:=^{self.width}}}'
        # title = self._trunc(title, self.width - 4)
        # print('', template.format(' ' + title + ' '), sep=os.linesep)
        print('')

    def h1(self, title):
        title = self._trunc(title, self.width)
        print('', title, '=' * len(title), sep=os.linesep)

    def h2(self, title):
        title = self._trunc(title, self.width)
        print('', title, '-' * len(title), sep=os.linesep)

    def describe(self, data):
        if data['type'] == 'hist_float':
            print('{}  count: {} '.format(data["name"], sum(data["counts"])),
                  'min: {:.8g} med: {:.8g} max: {:.8g}  '.format(data["min"], data["median"], data["max"]),
                  'std: {:.8g} var:{:.8g}'.format(data["std"], data["var"])
                  )

        elif data['type'] == 'hist_int':
            print('{}  count: {} '.format(data["name"], sum(data["counts"])),
                  'min: {:d} med: {:d} max: {:d}  '.format(int(data["min"]), int(data["median"]), int(data["max"]))
                  )

        elif data['type'] == 'hist_str':
            print('{} count: {} unique: {}'.format(data["name"], data["total"], data["unique"]))

        else:
            pass

    def hist(self, data: dict, width_hist=40):
        if data['type'] == 'hist_float':
            bin_edges = data['edges']

            ratio = width_hist / max(data['counts'])
            width_count = len(str(max(data['counts'])))

            for count, lower, upper in zip(data['counts'], bin_edges[:-1], bin_edges[1:]):
                scale = int(ratio * count)
                self.print('{:<{}} {:>{}d} [{:.4e}, {:.4e})'.format(
                    "▉" * scale, width_hist,
                    count, width_count,
                    lower, upper))

        elif data['type'] == 'hist_int':
            bin_edges = data['edges']

            ratio = width_hist / max(data['counts'])
            width_count = len(str(max(data['counts'])))

            for count, lower, upper in zip(data['counts'], bin_edges[:-1], bin_edges[1:]):
                scale = int(ratio * count)
                self.print('{:<{}} {:>{}d} [{:d}, {:d})'.format(
                    "▉" * scale, width_hist,
                    count, width_count,
                    np.ceil(lower).astype(int), np.floor(upper).astype(int)))

        elif data['type'] == 'hist_str':
            remain = data['total'] - sum(data['counts'])
            if remain > 0:
                data['counts'] = (*data['counts'], remain)
                data['labels'] = (*data['labels'], '...')

            width_count = len(str(max(data['counts'])))
            ratio = width_hist / max(data['counts'])
            for label, count in zip(data['labels'], data['counts']):
                scale = int(ratio * count)
                self.print('{:<{}} {:>{}d} {}'.format("▉" * scale, width_hist, count, width_count, label))

        elif data['type'] == 'hist_labels':

            width_count = len(str(max(data['counts'])))
            ratio = width_hist / max(data['counts'])
            for label, count in zip(data['labels'], data['counts']):
                scale = int(ratio * count)
                self.print('{:<{}} {:>{}d} {}'.format("▉" * scale, width_hist, count, width_count, label))

        else:
            pass

    @staticmethod
    def print(*args, **kwargs):
        print(*args, **kwargs)

    def _trunc(self, text, width=None):
        width = width if width else self.width
        return text if len(text) < width else text[:width - 3] + '...'
