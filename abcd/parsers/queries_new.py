# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar
# This program is distributed under the MIT License, see LICENSE.md.

import logging

from lark import Lark, Transformer
from lark.exceptions import LarkError

logger = logging.getLogger(__name__)


# TODO: Reversed operator in the grammar (value op prop VS prop op value VS IN)


class DebugTransformer(Transformer):  # pragma: no cover
    def __init__(self):
        super().__init__()

    def __default__(self, data, children, meta):
        print("Node: ", data, children)
        return data


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    with open("query_new.lark") as file:
        parser = Lark(file.read(), start="expression")

    transformer = DebugTransformer()

    queries = (
        " ",
        "single",
        "not single",
        "operator_gt > 23 ",
        "operator_gt > -2.31e-5 ",
        'string = "some string"',
        'regexp ~ ".*H"',
        "aa & not bb",
        "aa & bb > 23.54 | cc & dd",
        "aa and bb > 22 and cc > 33 and dd > 44 ",
        "((aa and bb > 22) and cc > 33) and dd > 44 ",
        "(aa and bb > 22) and (cc > 33 and dd > 44) ",
        "(aa and bb > 22 and cc > 33 and dd > 44) ",
        "aa and bb > 23.54 or 22 in cc and dd",
        "aa & bb > 23.54 | (22 in cc & dd)",
        "aa and bb > 23.54 or (22 in cc and dd)",
        "aa and not (bb > 23.54 or (22 in cc and dd))",
        "expression = (bb/3-1)*cc",
        "energy/n_atoms > 3",
        "1=3",
        "all(aa) > 3",
        "any(aa) > 3",
        "aa = False",
        "aa = [True, True, True]",
    )

    for query in queries:
        print(query)

        try:
            tree = parser.parse(query)
            # print(tree)
            # print(tree.pretty())
            print(transformer.transform(tree))

        except LarkError as err:
            raise NotImplementedError from err
