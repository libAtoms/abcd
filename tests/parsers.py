import unittest


class ParsingExtras(unittest.TestCase):
    """Testing extra argument parser """

    def setUp(self):
        from abcd.parsers.extras import parser
        self.parser = parser

    def test_empty(self):
        """Empty input
        """
        s = ('', ' ')
        out = self.parser.parse(' '.join(s))
        self.assertEqual(out, {})

    def test_single(self):
        """Single keys"""
        s = (
            ('flag', {'flag': True}),
            ('hyphen-ated', {'hyphen-ated': True}),
            (' start_with_a_separator', {'start_with_a_separator': True}),
            ('multiple_separators      ', {'multiple_separators': True}),
        )

        for string, expected in s:
            with self.subTest(string=string):
                out = self.parser.parse(string)
                self.assertEqual(out, expected)

    def test_key_value_pairs(self):
        """Key value pairs"""
        s = (
            ('string="astring"', {'string': 'astring'}),
            ('string_sapces = "astring"', {'string_sapces': 'astring'}),
        )

        for string, expected in s:
            with self.subTest(string=string):
                out = self.parser.parse(string)
                self.assertEqual(out, expected)

    def test_string(self):
        """Strings"""
        s = (
            ('quoted_string="quoted value"', {'quoted_string': 'quoted value'}),
            (r'quoted_string_escaped="esc\"aped"', {'quoted_string_escaped': 'esc"aped'}),
        )

        for string, expected in s:
            with self.subTest(string=string):
                out = self.parser.parse(string)
                self.assertEqual(out, expected)

    def test_boolean(self):
        """Boolean type"""
        s = (
            ('true_value', {'true_value': True}),
            ('true_value_long = true', {'true_value_long': True}),
            ('false_value = F', {'false_value': False}),
        )

        for string, expected in s:
            with self.subTest(string=string):
                out = self.parser.parse(string)
                self.assertEqual(out, expected)

    def test_numbers(self):
        """Numbers"""
        s = (
            ('integer=22', {'integer': 22}),
            ('floating=1.1', {'floating': 1.1}),
            ('scientific_float=1.2e7', {'scientific_float': 1.2e7}),
            ('scientific_float_2=5e-6', {'scientific_float_2': 5e-6}),
        )

        for string, expected in s:
            with self.subTest(string=string):
                out = self.parser.parse(string)
                self.assertEqual(out, expected)

    def test_arrays(self):
        """Arrays"""
        s = (
            ('int_array={1 2 3}', {'int_array': [1, 2, 3]}),
            ('array_nested=[[1,2],[3,4]]', {'array_nested': [[1, 2], [3, 4]]}),  # gets flattented if not 3x3
            ('array_many_other_quotes=({[4 8 12]})', {'array_many_other_quotes': [[[4, 8, 12]]]}),
            ('array_boolean={T F T F}', {'array_boolean': [True, False, True, False]}),
            ('array_bool_commas=[T, T, F, T]', {'array_bool_commas': [True, True, False, True]}),
        )

        for string, expected in s:
            with self.subTest(string=string):
                out = self.parser.parse(string)
                self.assertEqual(out, expected)

    def test_composite(self):
        """Composite"""
        s = (
            (' ', {}),
            ('int_array={1 2 3}', {'int_array': [1, 2, 3]}),
            ('array_nested=[[1,2],[3,4]]', {'array_nested': [[1, 2], [3, 4]]}),
            ('floating=1.1', {'floating': 1.1}),
            ('string="astring"', {'string': 'astring'}),
            ('hyphen-ated', {'hyphen-ated': True}),
            (' start_with_a_separator', {'start_with_a_separator': True}),
            ('scientific_float=1.2e7', {'scientific_float': 1.2e7}),
            ('array_many_other_quotes=({[4 8 12]})', {'array_many_other_quotes': [[[4, 8, 12]]]}),
            ('array_boolean={T F T F}', {'array_boolean': [True, False, True, False]}),
            ('array_bool_commas=[T, T, F, T]', {'array_bool_commas': [True, True, False, True]}),
            ('trailing', {'trailing': True}),
        )

        composite_string = ' '.join(string for string, _ in s)

        composite_expected = {}
        for _, expected in s:
            composite_expected.update(expected)

        out = self.parser.parse(composite_string)
        self.assertEqual(out, composite_expected)

    @unittest.skip("known issues / future features ")
    def test_missing(self):
        s = (
            '2body=33.3',
            'Properties=species:S:1:pos:R:3',
            'double_equals=abc=xyz',

            u'\xfcnicode_key=val\xfce',
            u'unquoted_special_value=a_to_Z_$%%^&*\xfc\u2615',
            u'quoted_string_unicode="a_to_Z_$%%^&*\xfc\u2615"',

            'float_array="3.3 4.4"',
            'scientific_float_array="1.2 2.2e3 4e1 3.3e-1 2e-2"',
            'a3x3_array="1 4 7 2 5 8 3 6 9" '  # fortran ordering
            'Lattice="  4.3  0.0 0.0 0.0  3.3 0.0 0.0 0.0  7.0 " '  # spaces in array
            'comma_separated="7, 4, -1"',
            'array_boolean_2=" T, F, T " '  # leading spaces

            'not_array="1.2 3.4 text"',
            'not_bool_array=[T F S]',
        )
        pass


class ParsingQueries(unittest.TestCase):

    def setUp(self):
        from abcd.parsers.queries import parser
        self.parser = parser

    def test_operators(self):
        """Operators"""
        s = (
            ('single', {}),
            ('operator_gt > 23 ', {}),
            ('string = "some string"', {}),
            ('regexp ~ ".*H"', {}),
        )
        for string, expected in s:
            with self.subTest(string=string):
                self.parser.parse(string)

    def test_combination(self):
        """Combinations"""
        s = (
            ('aa & not bb', {}),
            ('aa & bb > 23.54 | cc & dd', {}),
            ('aa and bb > 23 and bb > 23 and bb > 23 ', {}),
            ('aa and bb > 23.54 or 22 in cc and dd', {}),
        )
        for string, expected in s:
            with self.subTest(string=string):
                self.parser.parse(string)

    def test_expressions(self):
        """Complex expressions and functions"""
        s = (
            ('expression = (bb/3-1)*cc', {}),
            ('energy/n_atoms > 3', {}),
            ('all(aa) > 3', {}),
            ('any(aa) > 3', {}),
        )
        for string, expected in s:
            with self.subTest(string=string):
                self.parser.parse(string)

    @unittest.skip("known issues / future features ")
    def test_missing(self):
        s = (
            ('1=3', {}),
            ('aa = [True True True]', {}),
            ('aa & bb > 23.54 | (22 in cc & dd)', {}),
            ('aa and bb > 23.54 or (22 in cc and dd)', {}),
            ('aa and (bb > 23.54 or (22 in cc and dd))', {}),
        )


if __name__ == '__main__':
    unittest.main()
