from abcd.query.parser import QueryParser


class MongoQuery(object):

    def __init__(self):
        pass

    def visit(self, syntax_tree):
        op, *args = syntax_tree
        fun = self.__getattribute__(f'visit_{op.lower()}')
        return fun(*args)

    def visit_name(self, fields):
        return {fields: {'$exist': True}}

    def visit_and(self, *args):
        return {'$and': [self.visit(arg) for arg in args]}
        # TODO recursively combining all the and statements
        # out = {}
        # for arg in args:
        #     a = self.visit(arg)
        #
        #     out.update(**a)
        # return out

    def visit_or(self, *args):
        return {'$or': [self.visit(arg) for arg in args]}

    def visit_eq(self, field, value):
        return {field[1]: value[1]}

    def visit_gt(self, field, value):
        return {field[1]: {'$gt': value[1]}}

    def visit_gte(self, field, value):
        return {field[1]: {'$gte': value[1]}}

    def visit_lt(self, field, value):
        return {field[1]: {'$lt': value[1]}}

    def visit_lte(self, field, value):
        return {field[1]: {'$lte': value[1]}}

    def visit_in(self, field, *values):
        return {field[1]: {'$in': [value[1] for value in values]}}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __call__(self, ast):
        out = {}

        if ast is None:
            return out

        return self.visit(ast)


if __name__ == '__main__':
    parser = QueryParser()

    queries = (
        '',
        'aa & bb > 23 ',
        'aa & bb > 23.54 | (22 in cc & dd)',
        'aa and a=2 and bb > 23.54 or (22 in cc and dd)',
        'aa and (bb > 23.54 or (22 in cc and dd))',
        'variable = "some string"',
    )

    # for query in queries:
    #     print(query)
    #     print(parser.parse(query))

    n = 3
    with MongoQuery() as model:
        print(queries[n])
        print(parser.parse(queries[n]))
        print(model(parser.parse(queries[n])))
