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
                self.print('{:<{}} {:>{}d} {:.2f} - {:.2f}'.format(
                    "▉" * scale, width_hist,
                    count, width_count,
                    lower, upper))

        elif data['type'] == 'hist_int':
            bin_edges = data['edges']

            ratio = width_hist / max(data['counts'])
            width_count = len(str(max(data['counts'])))

            for count, lower, upper in zip(data['counts'], bin_edges[:-1], bin_edges[1:]):
                scale = int(ratio * count)
                self.print('{:<{}} {:>{}d} {:d} - {:d}'.format(
                    "▉" * scale, width_hist,
                    count, width_count,
                    int(lower), int(upper)))

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
