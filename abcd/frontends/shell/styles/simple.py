import os
import logging
import numpy as np

from pprint import pprint
from tabulate import tabulate
from collections import Counter

from .abstract import Style

logger = logging.getLogger(__name__)


class SimpleStyle(Style):
    partialBlocks = ["▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"]  # char=pb

    def __init__(self, width=80):
        self.width = width

    def title(self, title):
        template = f'{{:=^{self.width}}}'
        title = self._trunc(title, self.width - 4)
        print('', template.format(' ' + title + ' '), sep=os.linesep)

    def h1(self, title):
        title = self._trunc(title, self.width)
        print('', title, '=' * len(title), sep=os.linesep)

    def h2(self, title):
        title = self._trunc(title, self.width)
        print('', title, '-' * len(title), sep=os.linesep)

    def describe(self, data):
        if data['type'] == 'hist_float':

            print(
                f'{data["name"]}  count: {sum(data["counts"])} '
                f'min: {data["min"]:.8g} med: {data["median"]:.8g} max: {data["max"]:.8g}  '
                f'std: {data["std"]:.8g} var:{data["var"]:.8g}'
            )

        elif data['type'] == 'hist_str':

            print(f'{data["name"]} count: {data["total"]} unique: {data["unique"]}')

        else:
            pass

    def hist(self, data: dict, width_hist=16):

        if data['type'] == 'hist_float':

            ratio = width_hist / max(data['counts'])
            width_count = len(str(max(data['counts'])))
            for count in data['counts']:
                scale = int(ratio * count)
                self.print(f'{"▉" * scale:<{width_hist}} {count:>{width_count}d}')

        elif data['type'] == 'hist_str':

            width_count = len(str(max(data['counts'])))
            ratio = width_hist / max(data['counts'])
            for label, count in zip(data['labels'], data['counts']):
                scale = int(ratio * count)
                self.print(f'{"▉" * scale:<{width_hist}} {count:>{width_count}d} {label}')

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
