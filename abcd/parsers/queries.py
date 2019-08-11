import logging
from lark import Lark, Transformer, v_args
from lark.exceptions import LarkError

logger = logging.getLogger(__name__)

# https://github.com/lark-parser/lark/blob/master/examples/calc.py
# https://github.com/Materials-Consortia/optimade-python-tools/tree/master/optimade/grammar

grammar = r"""
    start: expression |
    //expression: [ NOT ] NAME  [ OPERATOR value ] [ conjunction expression ]
    expression: [ NOT ] ( NAME | sum ) [ OPERATOR ( value | sum ) ] [ conjunction expression ] 
              | [ NOT ] ( value | sum ) [ REVERSED ( NAME | sum ) ] [ conjunction expression ]

    OPERATOR: "=" | "!=" | ">" | ">=" | "<" | "<="| "~" 
    REVERSED: "in"

    conjunction: AND | OR 
    AND: "&" | "and"
    OR: "|" | "or"
    NOT: "!" | "not"

    value: NAME | SIGNED_FLOAT | SIGNED_INT | ESCAPED_STRING 
    
    // arithmetic
    ?sum: product
        | func
        | sum "+" product   -> add
        | sum "-" product   -> sub
    ?product: atom
        | product "*" atom  -> mul
        | product "/" atom  -> div
    ?atom: NUMBER           -> number
         | "-" atom         -> neg
         | NAME             -> var
         | "(" sum ")"

    // functions
    func: "all" "(" NAME ")" | "any" "(" NAME ")"
         
    %import common.CNAME -> NAME
    %import common.SIGNED_FLOAT
    %import common.SIGNED_INT
    %import common.ESCAPED_STRING
    %import common.WS_INLINE
    %import common.NUMBER

    %ignore WS_INLINE
"""


class TreeToAST(Transformer):
    pass


# parser = Lark(grammar, parser='lalr', lexer='contextual', transformer=TreeToAST(), debug=False)
parser = Lark(grammar, debug=False)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    queries = (
        'single'
        'operator_gt > 23 ',
        'string = "some string"',
        'regexp ~ ".*H"',
        'aa & not bb',
        'aa & bb > 23.54 | cc & dd',
        'aa and bb > 23 and bb > 23 and bb > 23 ',
        'aa and bb > 23.54 or 22 in cc and dd',
        'expression = (bb/3-1)*cc',
        'energy/n_atoms > 3',
        '1=3',
        'all(aa) > 3',
        'any(aa) > 3',
        # 'aa = [True True True]',
        # 'aa & bb > 23.54 | (22 in cc & dd)',
        # 'aa and bb > 23.54 or (22 in cc and dd)',
        # 'aa and (bb > 23.54 or (22 in cc and dd))',
    )

    for query in queries:
        print(query)
        print(parser.parse(query).pretty())

    # if '':
    print(parser.parse(''))
