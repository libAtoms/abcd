import logging
from lark import Lark, Transformer, v_args
from lark.exceptions import LarkError

logger = logging.getLogger(__name__)

# https://github.com/lark-parser/lark/blob/master/examples/calc.py
# https://github.com/Materials-Consortia/optimade-python-tools/tree/master/optimade/grammar

grammar = r"""
    ?start: expression | 
    ?expression: (negation | single_expression | reversed_expression | grouped_expression) [ relation expression ] 

    single_expression:   ( NAME | function ) [ operators value ]
    reversed_expression: value reversed_operators NAME  
   
    grouped_expression:   "(" expression ")"
    negation: NOT (single_expression | reversed_expression | grouped_expression)
    
    ?operators: EQUAL
             | NOTEQUAL
             | GT
             | GTE
             | LT
             | LTE
             | REGEXP

    ?reversed_operators: IN
        
    EQUAL:    "=" 
    NOTEQUAL: "!=" 
    GT:       ">" 
    GTE:      ">=" 
    LT:       "<" 
    LTE:      "<="
    REGEXP:   "~"     
    IN:       "in"

    ?relation: AND | OR 

    AND: "&" | "and"
    OR:  "|" | "or"
    NOT: "!" | "not"

    value: NAME 
         | FLOAT  -> float
         | INT    -> int
         | STRING -> string
         | "True"             -> true
         | "False"            -> false
         | array

    array: "[" value ([","] value)* "]"
    
    function: all | any
    
    all: "all" "(" NAME ")" 
    any: "any" "(" NAME ")"

     
    %import common.CNAME          -> NAME
    %import common.SIGNED_FLOAT   -> FLOAT
    %import common.SIGNED_INT     -> INT
    %import common.ESCAPED_STRING -> STRING
    %import common.WS_INLINE

    %ignore WS_INLINE
"""


class TreeToAST(Transformer):
    def start(self, _):
        return None

    int = v_args(inline=True)(int)
    float = v_args(inline=True)(float)
    array = list

    true = lambda self, _: True
    false = lambda self, _: False

    @v_args(inline=True)
    def string(self, s):
        return s[1:-1].replace('\\"', '"')

    def single_expression(self, items):

        if len(items) == 1:
            return {str(items[0]): {'EXISTS': True}}

        if len(items) != 3:
            raise NotImplementedError

        return {str(items[0]): {items[1].type: items[2]}}

    def reversed_expression(self, items):

        if len(items) != 3:
            raise NotImplementedError

        return {str(items[2]): {items[1].type: items[0]}}

    def grouped_expression(self, items):
        return items[0]

    def negation(self, s):
        return {k: {s[0].type: v} for k, v in s[1].items()}

    def expression(self, items):

        if len(items) != 3:
            raise NotImplementedError

        return {items[1].type: [items[0], items[2]]}


class Parser:
    parser = Lark(grammar)

    # parser = Lark(grammar, parser='lalr', lexer='contextual', transformer=TreeToAST(), debug=False)

    def parse(self, string):
        return self.parser.parse(string)


parser = Parser()
transformer = TreeToAST()

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    queries = (
        ' ',
        'single',
        'not single',
        'operator_gt > 23 ',
        'operator_gt > -2.31e-5 ',
        'string = "some string"',
        'regexp ~ ".*H"',
        'aa & not bb',
        'aa & bb > 23.54 | cc & dd',
        'aa and bb > 23 and bb > 23 and bb > 23 ',
        'aa and bb > 23.54 or 22 in cc and dd',
        'aa & bb > 23.54 | (22 in cc & dd)',
        'aa and bb > 23.54 or (22 in cc and dd)',
        'aa and not (bb > 23.54 or (22 in cc and dd))',
        # # 'expression = (bb/3-1)*cc',
        # # 'energy/n_atoms > 3',
        # '1=3',
        # 'all(aa) > 3',
        # 'any(aa) > 3',
        # 'aa = False',
        # 'aa = [True True True]',
    )

    for query in queries:
        # print(query)

        # print(parser.parse(query).pretty())
        try:
            tree = parser.parse(query)
            # print('=> tree:', tree)
            print('==> ast:', transformer.transform(tree))
        except LarkError:
            raise NotImplementedError
