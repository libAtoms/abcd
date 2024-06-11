from collections import Counter
from datetime import datetime
import logging
from operator import itemgetter

import numpy as np

logger = logging.getLogger(__name__)


def histogram(name, data, **kwargs):
    if not data:
        return None

    if isinstance(data, list):
        ptype = type(data[0])

        if not all(isinstance(x, ptype) for x in data):
            print("Mixed type error of the %s property!", name)
            return None

        if ptype == float:
            bins = kwargs.get("bins", 10)
            return _hist_float(name, data, bins)

        if ptype == int:
            bins = kwargs.get("bins", 10)
            return _hist_int(name, data, bins)

        if ptype == str:
            return _hist_str(name, data, **kwargs)

        if ptype == datetime:
            bins = kwargs.get("bins", 10)
            return _hist_date(name, data, bins)

        print(
            "%s: Histogram for list of %s types are not supported!", name, type(data[0])
        )
        logger.info(
            "%s: Histogram for list of %s types are not supported!", name, type(data[0])
        )

    logger.info("%s: Histogram for %s types are not supported!", name, type(data))
    return None


def _hist_float(name, data, bins=10):
    data = np.array(data)
    hist, bin_edges = np.histogram(data, bins=bins)

    return {
        "type": "hist_float",
        "name": name,
        "bins": bins,
        "edges": bin_edges,
        "counts": hist,
        "min": data.min(),
        "max": data.max(),
        "median": data.mean(),
        "std": data.std(),
        "var": data.var(),
    }


def _hist_date(name, data, bins=10):
    hist_data = np.array([t.timestamp() for t in data])
    hist, bin_edges = np.histogram(hist_data, bins=bins)

    fromtimestamp = datetime.fromtimestamp

    return {
        "type": "hist_date",
        "name": name,
        "bins": bins,
        "edges": [fromtimestamp(d) for d in bin_edges],
        "counts": hist,
        "min": fromtimestamp(hist_data.min()),
        "max": fromtimestamp(hist_data.max()),
        "median": fromtimestamp(hist_data.mean()),
        "std": fromtimestamp(hist_data.std()),
        "var": fromtimestamp(hist_data.var()),
    }


def _hist_int(name, data, bins=10):
    data = np.array(data)
    delta = max(data) - min(data) + 1

    bins = min(bins, delta)

    hist, bin_edges = np.histogram(data, bins=bins)

    return {
        "type": "hist_int",
        "name": name,
        "bins": bins,
        "edges": bin_edges,
        "counts": hist,
        "min": data.min(),
        "max": data.max(),
        "median": data.mean(),
        "std": data.std(),
        "var": data.var(),
    }


def _hist_str(name, data, bins=10, truncate=20):
    n_unique = len(set(data))

    if truncate:
        # data = (item[:truncate] for item in data)
        data = (
            item[:truncate] + "..." if len(item) > truncate else item for item in data
        )

    data = Counter(data)

    if bins:
        labels, counts = zip(*sorted(data.items(), key=itemgetter(1, 0), reverse=True))
    else:
        labels, counts = zip(*data.items())

    return {
        "type": "hist_str",
        "name": name,
        "total": sum(data.values()),
        "unique": n_unique,
        "labels": labels[:bins],
        "counts": counts[:bins],
    }
