# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar
# This program is distributed under the MIT License, see LICENSE.md.

import sys
from lark import Lark, Transformer, v_args
from lark.exceptions import LarkError

grammar = r"""
    start: ( key | key_value )*

    key: NAME
    key_value: NAME "=" value

    NAME: ("_"|LETTER|DIGIT) ("_"|"-"|LETTER|DIGIT)*

    ?value: 
          | ESCAPED_STRING     -> string
          | integer
          | float
          | boolean
          | WORD
          | "null"             -> null
          
    arrays: "[" value (","? value)* "]"
          | "{" value (","? value)* "}"
          | "(" value (","? value)* ")"
    
    ?integer: INT -> int 
            | array_int
    
    ?float: FLOAT -> float 
          | array_float
    
    ?boolean: true | false | array_boolean

    true: "T" | "true"
    false: "F"| "false"

    array_int: "[" integer ([","] integer)* "]"
             | "{" integer ([","] integer)* "}"
             | "(" integer ([","] integer)* ")"
             
    
    array_float: "[" float ([","] float)* "]"
               | "{" float ([","] float)* "}"
               | "(" float ([","] float)* ")"

    array_boolean: "[" boolean ([","] boolean)* "]"
               | "{" boolean ([","] boolean)* "}"
               | "(" boolean ([","] boolean)* ")"
    
    
    INT.3: SIGNED_INT 
    FLOAT.4: SIGNED_FLOAT
          
    %import common.ESCAPED_STRING
    %import common.LETTER
    %import common.WORD
    %import common.DIGIT
    %import common.SIGNED_NUMBER
    %import common.SIGNED_INT
    %import common.SIGNED_FLOAT
    %import common.WS_INLINE
    %ignore WS_INLINE           // Disregard spaces in text
"""


class TreeToDict(Transformer):
    @v_args(inline=True)
    def string(self, s):
        return s[1:-1].replace('\\"', '"')

    arrays = list
    pair = tuple
    object = dict
    number = v_args(inline=True)(float)

    int = lambda self, s: int(s[0])
    float = lambda self, s: float(s[0])
    true = lambda self, _: True
    false = lambda self, _: False
    null = lambda self, _: None

    # key_value = tuple
    key_value = lambda self, s: (str(s[0]), s[1])
    key = lambda self, s: (str(s[0]), True)
    start = dict

    namee = print
    array_int = list
    array_float = list
    array_boolean = list


# parser = Lark(grammar, parser='lalr', lexer='contextual', debug=False)
parser = Lark(
    grammar, parser="lalr", lexer="contextual", transformer=TreeToDict(), debug=False
)

if __name__ == "__main__":
    test_string = " ".join(
        [
            " " "flag",  # start with a separator
            'quotedd_string="quoteddd value"',
            r'quotedddd_string_escaped="esc\"aped"',
            "false_value = F",
            "integer=22",
            "floating=1.1",
            "int_array={1 2 3}",
            "scientific_float=1.2e7",
            "scientific_float_2=5e-6",
            'scientific_float_array="1.2 2.2e3 4e1 3.3e-1 2e-2"',
            'not_array="1.2 3.4 text"',
            "array_nested=[[1,2],[3,4]] "  # gets flattented if not 3x3
            "array_many_other_quotes=({[4 8 12]})",
            "array_boolean={T F T F}",
            'array_boolean_2=" T, F, T " '  # leading spaces
            # 'not_bool_array=[T F S]',
            # # read and write
            # u'\xfcnicode_key=val\xfce',
            # # u'unquoted_special_value=a_to_Z_$%%^&*\xfc\u2615',
            # '2body=33.3',
            "hyphen-ated",
            # # # parse only
            'comma_separated="7, 4, -1"',
            "array_bool_commas=[T, T, F, T]",
            # # 'Properties=species:S:1:pos:R:3',
            # 'double_equals=abc=xyz',
            "multiple_separators      ",
            "trailing",
        ]
    )

    j = parser.parse(test_string)

    print(j)
    # TODO: arrays with mixed types
    # print(j.pretty())

    # try:
    #     parser.parse('double_equals=abc=xyz')
    # except LarkError:
    #     print('ERROR! ' * 10)

# #
# # This example shows how to write a basic JSON parser
# #
# # The code is short and clear, and outperforms every other parser (that's written in Python).
# # For an explanation, check out the JSON parser tutorial at /docs/json_tutorial.md
# #
#
# import sys
#
# from lark import Lark, Transformer, v_args
#
# json_grammar = r"""
#     ?start: value
#
#     ?value: object
#           | array
#           | string
#           | SIGNED_NUMBER      -> number
#           | "true"             -> true
#           | "false"            -> false
#           | "null"             -> null
#
#     array  : "[" [value ("," value)*] "]"
#     object : "{" [pair ("," pair)*] "}"
#     pair   : string ":" value
#
#     string : ESCAPED_STRING
#
#     %import common.ESCAPED_STRING
#     %import common.SIGNED_NUMBER
#     %import common.WS
#
#     %ignore WS
# """
#
#
# class TreeToJson(Transformer):
#     @v_args(inline=True)
#     def string(self, s):
#         return s[1:-1].replace('\\"', '"')
#
#     array = list
#     pair = tuple
#     object = dict
#     number = v_args(inline=True)(float)
#
#     null = lambda self, _: None
#     true = lambda self, _: True
#     false = lambda self, _: False
#
#
# # json_parser = Lark(json_grammar, parser='earley', lexer='standard')
# # def parse(x):
# #     return TreeToJson().transform(json_parser.parse(x))
#
# json_parser = Lark(json_grammar, parser='lalr', lexer='standard', transformer=TreeToJson())
# parse = json_parser.parse
#
#
# def test():
#     test_json = '''
#         {
#             "empty_object" : {},
#             "empty_array"  : [],
#             "booleans"     : { "YES" : true, "NO" : false },
#             "numbers"      : [ 0, 1, -2, 3.3, 4.4e5, 6.6e-7 ],
#             "strings"      : [ "This", [ "And" , "That", "And a \\"b" ] ],
#             "nothing"      : null
#         }
#     '''
#
#     j = parse(test_json)
#     print(j)
#     import json
#     assert j == json.loads(test_json)
#
#
# if __name__ == '__main__':
#     # test()
#     with open(sys.argv[1]) as f:
#         print(parse(f.read()))
