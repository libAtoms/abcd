from abcd.query.parser import QueryParser

if __name__ == '__main__':
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
        print(parser.parse(query))
