import logging
import ply.lex as lex
import ply.yacc as yacc

logger = logging.getLogger(__name__)

reserved = {
    'and': 'AND',
    'or': 'OR',
    'in': 'IN',
    'not_in': 'NIN',
}


# noinspection PySingleQuotedDocstring
class QueryLexer(object):
    # List of token names.   This is always required
    tokens = ('FLOAT', 'INT', 'NAME', 'STRING',
              'EQ', 'RE', 'LT', 'GT', 'LTE', 'GTE',
              'LPAREN', 'RPAREN') + tuple(reserved.values())

    # tokens = ('FLOAT', 'INT', 'NAME', 'STRING',
    #           'EQ', 'LT', 'GT', 'LTE', 'GTE',
    #           'LPAREN', 'RPAREN', 'BEGIN_ARRAY', 'END_ARRAY', 'VALUE_SEPARATOR') + tuple(reserved.values())

    # Tokens
    t_FLOAT = r'[-+]?[0-9]+(\.([0-9]+)?([eE][-+]?[0-9]+)?|[eE][-+]?[0-9]+)'
    t_INT = r'[-+]?[0-9]+'

    t_STRING = r'"[a-zA-Z0-9_ ]*"'

    # Special symbols/operators
    t_AND = r'&'
    t_OR = r'\|'

    t_EQ = r'=+'
    t_RE = r'~+'

    t_LT, t_LTE = r'<', r'<='
    t_GT, t_GTE = r'>', r'>='

    t_LPAREN, t_RPAREN = r'\(', r'\)'

    # t_BEGIN_ARRAY, t_END_ARRAY = r'\[', r'\]'
    # t_VALUE_SEPARATOR = r','

    # A string containing ignored characters (spaces and tabs)
    t_ignore = ' \t'

    def __init__(self, **kwargs):
        # Build the lexer
        self.lexer = lex.lex(module=self, **kwargs)

    @staticmethod
    def t_name(t):
        r'[a-zA-Z_\-.^\[\]][a-zA-Z0-9_\-.*^\[\]]*[$]?'
        t.type = reserved.get(t.value, 'NAME')  # Check for reserved words
        return t

    # Define a rule so we can track line numbers
    @staticmethod
    def t_newline(t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    # Error handling rule
    @staticmethod
    def t_error(t):
        logger.warning('Illegal character \'{}\''.format(t.value[0]))
        t.lexer.skip(1)

    def tokenize(self, data):
        self.lexer.input(data)
        return self.lexer


class QueryParser(QueryLexer):
    """
    Base class for a lexer/parser that has the rules defined as methods
    """
    precedence = (
        ('left', 'AND', 'OR', 'IN', 'NIN'),
        ('left', 'EQ', 'RE', 'LT', 'GT', 'LTE', 'GTE'),
    )

    def __init__(self, **kwargs):
        # Build the lexer and parser
        super().__init__(**kwargs)
        self.parser = yacc.yacc(module=self, debug=False, write_tables=False, **kwargs)

    @staticmethod
    def p_statements_single(p):
        """statements : statement"""
        p[0] = p[1]

    @staticmethod
    def p_statements_multiple(p):
        """statements : statements statement"""

        p[0] = ('AND',
                p[1][1:] if p[1][0] == 'AND' else p[1],
                p[2][1:] if p[2][0] == 'AND' else p[2])

    @staticmethod
    def p_statement_expr(p):
        """statement : expression"""
        p[0] = p[1]

    @staticmethod
    def p_expression_int(p):
        """expression : INT"""
        p[0] = ('NUMBER', int(p[1]))

    @staticmethod
    def p_expression_float(p):
        """expression : FLOAT"""
        p[0] = ('NUMBER', float(p[1]))

    @staticmethod
    def p_expression_string(p):
        """expression : STRING"""
        p[0] = ('STRING', p[1].strip('"'))

    @staticmethod
    def p_expression_name(p):
        """expression : NAME"""
        p[0] = ('NAME', p[1])

    # @staticmethod
    # def p_expression_name_exist(p):
    #     """expression : NAME expression"""
    #     p[0] = ('AND', ('EXSITS', p[1]), p[2])

    @staticmethod
    def p_expression_and(p):
        """expression : expression AND expression"""
        out = ('AND',)
        if p[1][0] == 'AND':
            out += p[1][1:]
        else:
            out += (p[1],)

        if p[3][0] == 'AND':
            out += p[3][1:]
        else:
            out += (p[3],)

        p[0] = out

    @staticmethod
    def p_expression_binop(p):
        """expression : expression OR expression
                      | expression IN expression
                      | expression NIN expression
                      | expression EQ expression
                      | expression RE expression
                      | expression LT expression
                      | expression GT expression
                      | expression LTE expression
                      | expression GTE expression
        """
        # p[0] = (p[2], p[1], p[3])
        p[0] = (p.slice[2].type, p[1], p[3])

    @staticmethod
    def p_expression_group(p):
        """expression : LPAREN expression RPAREN"""
        p[0] = p[2]

    #
    # # @staticmethod
    # # def p_number(p):
    # #     """number : FLOAT"""
    # #     p[0] = p[1]
    #
    # @staticmethod
    # def p_value(p):
    #     """value : array
    #              | FLOAT
    #              | INT
    #     """
    #     p[0] = p[1]
    #
    # @staticmethod
    # def p_values(p):
    #     """values : values value VALUE_SEPARATOR
    #               | values value"""
    #     if len(p) == 1:
    #         p[0] = list()
    #     else:
    #         p[1].append(p[2])
    #         p[0] = p[1]
    #
    # def p_array(self, p):
    #     """array : BEGIN_ARRAY values END_ARRAY"""
    #     p[0] = p[2]

    @staticmethod
    def p_error(p):
        if p:
            logger.error('Syntax error at \'{}\''.format(p.value))
        else:
            logger.error("Syntax error at EOF")

    def parse(self, s):
        if not s:
            return None
        elif isinstance(s, dict):
            return s
        elif isinstance(s, str):
            return self.parser.parse(s)
        elif isinstance(s, list):
            return self.parser.parse(' and '.join(s))
        else:
            raise NotImplementedError('Query string should be string or dictionary and not {}!'.format(type(query)))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    lexer = QueryLexer()
    parser = QueryParser()

    queries = (
        # 'aa = [True True True]',
        'formula~.*H',
        'aa bb > 23 ',
        'aa-bb > 23 ',
        'aa & bb > 23 & bb > 23 & bb > 23 ',
        'aa & bb > 23.54 | (22 in cc & dd)',
        'aa and bb > 23.54 or (22 in cc and dd)',
        'aa and (bb > 23.54 or (22 in cc and dd))',
        'variable = "some string"',
    )

    for query in queries:
        result = lexer.tokenize(query)
        print(query)
        print(' '.join('{} [{}]'.format(item.value, item.type) for item in result))
        print(parser.parse(query))

    print(parser.parse(['aa bb > 23', 'aa and (bb > 23.54 or (22 in cc and dd))']))
