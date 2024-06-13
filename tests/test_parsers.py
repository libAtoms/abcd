import pytest
from abcd.parsers.extras import parser as extras_parser
from abcd.parsers.queries import parser as queries_parser


class TestParsingExtras:
    @pytest.fixture
    def parser(self):
        return extras_parser

    @pytest.mark.parametrize(
        "string, expected",
        [
            ("flag", {"flag": True}),
            ("hyphen-ated", {"hyphen-ated": True}),
            (" start_with_a_separator", {"start_with_a_separator": True}),
            ("multiple_separators      ", {"multiple_separators": True}),
        ],
    )
    def test_single(self, parser, string, expected):
        """Single keys"""
        assert expected == parser.parse(string)

    @pytest.mark.parametrize(
        "string, expected",
        [
            ('string="astring"', {"string": "astring"}),
            ('string_sapces = "astring"', {"string_sapces": "astring"}),
        ],
    )
    def test_key_value_pairs(self, parser, string, expected):
        """Key value pairs"""
        assert expected == parser.parse(string)

    @pytest.mark.parametrize(
        "string, expected",
        [
            ('quoted_string="quoted value"', {"quoted_string": "quoted value"}),
            (
                r'quoted_string_escaped="esc\"aped"',
                {"quoted_string_escaped": 'esc"aped'},
            ),
        ],
    )
    def test_string(self, parser, string, expected):
        """Strings"""
        assert expected == parser.parse(string)

    @pytest.mark.parametrize(
        "string, expected",
        [
            ("true_value", {"true_value": True}),
            ("true_value_long = true", {"true_value_long": True}),
            ("false_value = F", {"false_value": False}),
            ("false_value_colon: F", {"false_value_colon": False}),
        ],
    )
    def test_boolean(self, parser, string, expected):
        """Boolean type"""
        assert expected == parser.parse(string)

    @pytest.mark.parametrize(
        "string, expected",
        [
            ("integer=22", {"integer": 22}),
            ("floating=1.1", {"floating": 1.1}),
            ("scientific_float=1.2e7", {"scientific_float": 1.2e7}),
            ("scientific_float_2=5e-6", {"scientific_float_2": 5e-6}),
            ("floating_colon: 3.14", {"floating_colon": 3.14}),
        ],
    )
    def test_numbers(self, parser, string, expected):
        """Numbers"""
        assert expected == parser.parse(string)

    @pytest.mark.parametrize(
        "string, expected",
        [
            ("int_array={1 2 3}", {"int_array": [1, 2, 3]}),
            ("array_nested=[[1,2],[3,4]]", {"array_nested": [[1, 2], [3, 4]]}),
            # gets flattented if not 3x3
            (
                "array_many_other_quotes=({[4 8 12]})",
                {"array_many_other_quotes": [[[4, 8, 12]]]},
            ),
            ("array_boolean={T F T F}", {"array_boolean": [True, False, True, False]}),
            (
                "array_bool_commas=[T, T, F, T]",
                {"array_bool_commas": [True, True, False, True]},
            ),
            ("int_array_colon: {4 2}", {"int_array_colon": [4, 2]}),
        ],
    )
    def test_arrays(self, parser, string, expected):
        """Arrays"""
        assert expected == parser.parse(string)

    def test_composite(self, parser):
        """Composite"""
        s = (
            (" ", {}),
            ("int_array={1 2 3}", {"int_array": [1, 2, 3]}),
            ("array_nested=[[1,2],[3,4]]", {"array_nested": [[1, 2], [3, 4]]}),
            ("floating=1.1", {"floating": 1.1}),
            ('string="astring"', {"string": "astring"}),
            ("hyphen-ated", {"hyphen-ated": True}),
            (" start_with_a_separator", {"start_with_a_separator": True}),
            ("scientific_float=1.2e7", {"scientific_float": 1.2e7}),
            (
                "array_many_other_quotes=({[4 8 12]})",
                {"array_many_other_quotes": [[[4, 8, 12]]]},
            ),
            ("array_boolean={T F T F}", {"array_boolean": [True, False, True, False]}),
            (
                "array_bool_commas=[T, T, F, T]",
                {"array_bool_commas": [True, True, False, True]},
            ),
            ("trailing", {"trailing": True}),
        )

        composite_string = " ".join(string for string, _ in s)

        composite_expected = {}
        for _, expected in s:
            composite_expected.update(expected)

        out = parser.parse(composite_string)
        assert out == composite_expected

    @pytest.mark.parametrize(
        "string, expected",
        [
            ('colon_string:"astring"', {"colon_string": "astring"}),
            ('colon_string_spaces : "astring"', {"colon_string_spaces": "astring"}),
        ],
    )
    def test_colon_key_value_pairs(self, parser, string, expected):
        """Key value pairs separated by colons"""
        assert expected == parser.parse(string)

    @pytest.mark.skip
    @pytest.mark.parametrize(
        "string",
        [
            "2body=33.3",
            "Properties=species:S:1:pos:R:3",
            "double_equals=abc=xyz",
            "\xfcnicode_key=val\xfce",
            "unquoted_special_value=a_to_Z_$%%^&*\xfc\u2615",
            'quoted_string_unicode="a_to_Z_$%%^&*\xfc\u2615"',
            'float_array="3.3 4.4"',
            'scientific_float_array="1.2 2.2e3 4e1 3.3e-1 2e-2"',
            'a3x3_array="1 4 7 2 5 8 3 6 9" '  # fortran ordering
            'Lattice="  4.3  0.0 0.0 0.0  3.3 0.0 0.0 0.0  7.0 " '  # spaces in array
            'comma_separated="7, 4, -1"',
            'array_boolean_2=" T, F, T " ' 'not_array="1.2 3.4 text"',  # leading spaces
            "not_bool_array=[T F S]",
        ],
    )
    def test_missing(self, string):
        ...


class TestParsingQueries:
    @pytest.fixture
    def parser(self):
        return queries_parser

    @pytest.mark.parametrize(
        "string, expected",
        [
            ("single", {}),
            ("operator_gt > 23 ", {}),
            ('string = "some string"', {}),
            ('regexp ~ ".*H"', {}),
        ],
    )
    def test_operators(self, parser, string, expected):
        """Operators"""
        assert parser.parse(string)

    @pytest.mark.parametrize(
        "string, expected",
        [
            ("aa && not bb", {}),
            ("aa && bb > 23.54 || cc && dd", {}),
            ("aa and bb > 23 and bb > 23 and bb > 23 ", {}),
            ("aa and bb > 23.54 or 22 in cc and dd", {}),
        ],
    )
    def test_combination(self, parser, string, expected):
        """Combinations"""
        assert parser.parse(string)

    @pytest.mark.skip
    @pytest.mark.parametrize(
        "case",
        [
            ("expression = (bb/3-1)*cc", {}),
            ("energy/n_atoms > 3", {}),
            ("all(aa) > 3", {}),
            ("any(aa) > 3", {}),
        ],
    )
    def test_expressions(self, case):
        ...

    @pytest.mark.skip("known issues / future features")
    @pytest.mark.parametrize(
        "case",
        [
            ("1=3", {}),
            ("aa = [True True True]", {}),
            ("aa && bb > 23.54 || (22 in cc && dd)", {}),
            ("aa and bb > 23.54 or (22 in cc and dd)", {}),
            ("aa and (bb > 23.54 or (22 in cc and dd))", {}),
        ],
    )
    def test_expressions(self, case):
        ...
