import logging
from lark import Lark, Transformer, v_args
from lark.exceptions import LarkError
from abcd.queryset import Query

logger = logging.getLogger(__name__)

# https://github.com/lark-parser/lark/blob/master/examples/calc.py
# https://github.com/Materials-Consortia/optimade-python-tools/tree/master/optimade/grammar

grammar = r"""
    ?statement:                                   -> empty
              | expression                        -> single_statement
              | statement relation_ops expression -> multi_statement
    
    ?expression: NAME                      -> single_expression 
               | NAME comparison_ops value -> operator_expression 
               | value reversed_ops NAME   -> reversed_expression 
               | "(" statement ")"         -> grouped_expression 
               | NOT statement             -> negation_expression
        
    ?comparison_ops: EQUAL
              | NOTEQUAL
              | GT
              | GTE
              | LT
              | LTE
              | RE

        ?relation_ops: AND 
                 | OR 

    AND: "&" | "and"
    OR:  "|" | "or"
    NOT: "!" | "not"

    ?reversed_ops: IN
        
    EQUAL:    "=" 
    NOTEQUAL: "!=" 
    GT:       ">" 
    GTE:      ">=" 
    LT:       "<" 
    LTE:      "<="
    RE:       "~"     
    IN:       "in"


    value: FLOAT    -> float
         | INT      -> int
         | STRING   -> string
         | "True"   -> true
         | "False"  -> false
         | array

    array: "[" value ([","] value)* "]"
    
    //function: all | any
    
    //all: "all" "(" NAME ")" 
    //any: "any" "(" NAME ")"

    //NAME : /(?!and\b)/ /(?!or\b)/ /(?!not\b)/ CNAME

    %import common.CNAME          -> NAME          
    %import common.SIGNED_FLOAT   -> FLOAT
    %import common.SIGNED_INT     -> INT
    %import common.ESCAPED_STRING -> STRING
    %import common.WS_INLINE

    %ignore WS_INLINE
"""


@v_args(inline=True)
class TreeTransformer(Transformer):
    empty = lambda self: ()

    @v_args(inline=False)
    def array(self, *items):
        return list(items)

    true = lambda _: ('VALUE', True)
    false = lambda _: ('VALUE', False)

    def float(self, number):
        return 'NUMBER', float(number)

    def int(self, number):
        return 'NUMBER', int(number)

    def string(self, s):
        return 'STRING', s[1:-1].replace('\\"', '"')

    def single_statement(self, expression):
        return expression

    def multi_statement(self, statement, operator, expression):
        return operator.type, statement, expression

    def single_expression(self, name):
        return 'NAME', str(name)

    def grouped_expression(self, statement):
        # return statement
        return 'GROUP', statement

    def operator_expression(self, name, operator, value):
        return operator.type, ('NAME', str(name)), value

    def reversed_expression(self, value, operator, name):
        return operator.type, ('NAME', str(name)), value

    def negation_expression(self, operator, expression):
        return operator.type, expression


class Parser:
    parser = Lark(grammar, start='statement')
    transformer = TreeTransformer()

    def parse(self, string):
        return self.parser.parse(string)

    def __call__(self, string):
        return self.transformer.transform(self.parse(string))


parser = Parser()

if __name__ == '__main__':
    # logging.basicConfig(level=logging.DEBUG)
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
        # 'aa bb > 22 cc > 33 dd > 44 ',
        'aa and bb > 22 and cc > 33 and dd > 44 ',
        '((aa and bb > 22) and cc > 33) and dd > 44 ',
        '(aa and bb > 22) and (cc > 33 and dd > 44) ',
        '(aa and bb > 22 and cc > 33 and dd > 44) ',
        'aa and bb > 23.54 or 22 in cc and dd',
        'aa & bb > 23.54 | (22 in cc & dd)',
        'aa and bb > 23.54 or (22 in cc and dd)',
        'aa and not (bb > 23.54 or (22 in cc and dd))',
        # 'expression = (bb/3-1)*cc',
        # 'energy/n_atoms > 3',
        # '1=3',
        # 'all(aa) > 3',
        # 'any(aa) > 3',
        # 'aa = False',
        # 'aa = [True True True]',
    )

    for query in queries:
        logger.info(query)

        # print(parser.parse(query).pretty())
        try:
            tree = parser.parse(query)
            logger.info('=> tree: {}'.format(tree))
            logger.info('==> ast: {}'.format(parser(query)))
        except LarkError:
            raise NotImplementedError
