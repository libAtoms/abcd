import logging
from sly import Lexer, Parser

logger = logging.getLogger(__name__)

reserved = {
    'and': 'AND',
    'or': 'OR',
    'in': 'IN',
    'not_in': 'NIN',
}


class QueryLexer(Lexer):
    tokens = {'FLOAT', 'INT', 'NAME', 'STRING',
              'EQ', 'LT', 'GT', 'LTE', 'GTE',
              'LPAREN', 'RPAREN'}.union(reserved.values())
    ignore = ' \t'

    # Tokens
    FLOAT = r'[-+]?[0-9]+(\.([0-9]+)?([eE][-+]?[0-9]+)?|[eE][-+]?[0-9]+)'
    INT = r'[-+]?[0-9]+'

    NAME = r'[a-zA-Z_][a-zA-Z0-9_]*'
    STRING = r'"[a-zA-Z0-9_ ]*"'

    # Special symbols/operators
    AND = r'&'
    OR = r'\|'

    EQ = r'='
    LT, LTE = r'<', r'<='
    GT, GTE = r'>', r'>='

    LPAREN, RPAREN = r'\(', r'\)'

    # Ignored pattern
    ignore_newline = r'\n+'

    def NAME(self, t):
        t.type = reserved.get(t.value, 'NAME')  # Check for reserved words
        return t

    # Extra action for newlines
    def ignore_newline(self, t):
        self.lineno += t.value.count('\n')
        raise ValueError('Multiline expressions are not supported yet!')

    def error(self, t):
        # raise ValueError(f'Illegal character \'{t.value[0]}\'')
        print("Illegal character '%s'" % t.value[0])
        self.index += 1


class QueryParser(Parser):
    tokens = QueryLexer.tokens

    precedence = (
        ('left', 'AND', 'OR', 'IN', 'NIN'),
        ('left', 'EQ', 'LT', 'GT', 'LTE', 'GTE'),
    )

    # @_('NAME EQ expr')
    # def statement(self, p):
    #     self.names[p.NAME] = p.expr

    @_('expr')
    def exprs(self, p):
        return p.expr

    @_('exprs expr')
    def exprs(self, p):
        return p.exprs + p.expr

    # @_('statement')
    # def statements(self, p):
    #     return p.statement
    #
    # @_('statements statement')
    # def statements(self, p):
    #     return p.statements + [p.statement]
    #
    # @_('expr')
    # def statement(self, p):
    #     return p.expr

    @_('FLOAT')
    def expr(self, p):
        return 'NUMBER', float(p.FLOAT)

    @_('INT')
    def expr(self, p):
        return 'NUMBER', int(p.INT)

    @_('STRING')
    def expr(self, p):
        return 'STRING', p.STRING.strip('"')

    @_('NAME')
    def expr(self, p):
        return 'NAME', p.NAME

    @_('expr AND expr',
       'expr OR expr',
       'expr IN expr',
       'expr NIN expr',
       'expr EQ expr',
       'expr LT expr',
       'expr GT expr',
       'expr LTE expr',
       'expr GTE expr')
    def expr(self, p):
        return p._slice[1].type, p.expr0, p.expr1

    @_('LPAREN expr RPAREN')
    def expr(self, p):
        return p.expr

    def error(self, t):
        # raise ValueError(f'Illegal character \'{t.value[0]}\'')
        print("Illegal character '%s'" % t.value[0])
        self.index += 1


if __name__ == '__main__':
    lexer = QueryLexer()
    parser = QueryParser()

    queries = (
        'aa & bb > 23 ',
        'aa & bb > 23.54 | (22 in cc & dd)',
        'aa and bb > 23.54 or (22 in cc and dd)',
        'aa and (bb > 23.54 or (22 in cc and dd))',
        'variable = "some string"',
    )
    for query in queries:
        print(query)
        print(' '.join(f'{item.value} [{item.type}]' for item in lexer.tokenize(query)))
        print(parser.parse(lexer.tokenize(query)))
