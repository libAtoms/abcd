import os
import shutil
import logging
import numpy as np

from pprint import pprint
from collections import Counter

from tabulate import tabulate

logger = logging.getLogger(__name__)

from .abstract import Style


class FancyStyle(Style):
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

    def describe(self, data):
        if data['type'] == 'hist_float':

            self.title(f'Histogram: {data["name"]}')
            self.h1(f'Property: {data["name"]}')

            self.h2('Summary')
            self.table(
                ((f'{sum(data["counts"])}', f'{data["min"]:.8g}', f'{data["max"]:.8g}', f'{data["median"]:.8g}',
                  f'{data["std"]:.8g}',
                  f'{data["var"]:.8g}'),),
                headers=('total', 'min', 'max', 'mean', 'std', 'var'), disable_numparse=True, tablefmt="grid")

        elif data['type'] == 'hist_str':

            self.title(f'Histogram: {data["name"]}')
            self.h1(f'Property: {data["name"]}')

            self.h2('Summary')
            self.table(
                ((f'{data["total"]}', f'{data["unique"]}', f'{min(data["counts"])}', f'{max(data["counts"])}',
                  f'{data["total"] / data["unique"]:.8g}'),),
                headers=('total', 'unique', 'min', 'max', 'mean'), disable_numparse=True, tablefmt="grid")
        else:
            pass

    def hist(self, data, **kwargs):
        if data['type'] == 'hist_float':

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

            hist, bin_edges = data['counts'], data['edges']
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

        elif data['type'] == 'hist_str':

            self.h2('Histogram')

            def get_width_hack(keys, values):
                # generating a table just to measure its width
                return len(tabulate(
                    list([item, f'{count}'] for item, count in zip(data['labels'], data['counts'])),
                    headers=('Name', 'Count'),
                    colalign=("left", "right"),
                    disable_numparse=True,
                    tablefmt="psql"
                ).split('\n', 1)[0])

            table_width = get_width_hack(data['labels'], data['counts'])

            if table_width > self.width:
                self.print('Too small terminal!')
            else:
                max_width = self.width - table_width - 3

                ratio = max_width / max(data['counts'])
                scales = (int(ratio * value) for value in data['counts'])

                self.table(
                    list([item, f'{count}', '▉' * scale] for item, count, scale in
                         zip(data['labels'], data['counts'], scales)),
                    headers=('Name', 'Count', 'Histogram'),
                    colalign=("left", "right", "left"),
                    disable_numparse=True,
                    tablefmt="psql"
                )
        elif data['type'] == 'hist_labels':

            def get_width_hack(keys, values):
                # generating a table just to measure its width
                return len(tabulate(
                    list([item, f'{count}'] for item, count in zip(data['labels'], data['counts'])),
                    headers=('Name', 'Count'),
                    colalign=("left", "right"),
                    disable_numparse=True,
                    tablefmt="psql"
                ).split('\n', 1)[0])

            table_width = get_width_hack(data['labels'], data['counts'])

            if table_width > self.width:
                self.print('Too small terminal!')
            else:
                max_width = self.width - table_width - 3

                ratio = max_width / max(data['counts'])
                scales = (int(ratio * value) for value in data['counts'])

                self.table(
                    list([item, f'{count}', '▉' * scale] for item, count, scale in
                         zip(data['labels'], data['counts'], scales)),
                    headers=('Name', 'Count', 'Histogram'),
                    colalign=("left", "right", "left"),
                    disable_numparse=True,
                    tablefmt="psql"
                )
        else:
            pass

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
