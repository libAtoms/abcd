import re
import numpy as np

UNPROCESSED_KEYS = ['uid']


def key_val_str_to_dict(string, sep=None):
    """
    Parse an xyz properties string in a key=value and return a dict with
    various values parsed to native types.

    Accepts brackets or quotes to delimit values. Parses integers, floats
    booleans and arrays thereof. Arrays with 9 values are converted to 3x3
    arrays with Fortran ordering.

    If sep is None, string will split on whitespace, otherwise will split
    key value pairs with the given separator.

    """
    # store the closing delimiters to match opening ones
    delimiters = {
        "'": "'",
        '"': '"',
        '(': ')',
        '{': '}',
        '[': ']',
    }

    # Make pairs and process afterwards
    kv_pairs = [
        [[]]]  # List of characters for each entry, add a new list for new value
    delimiter_stack = []  # push and pop closing delimiters
    escaped = False  # add escaped sequences verbatim

    # parse character-by-character unless someone can do nested brackets
    # and escape sequences in a regex
    for char in string.strip():
        if escaped:  # bypass everything if escaped
            kv_pairs[-1][-1].extend(['\\', char])
            escaped = False
        elif delimiter_stack:  # inside brackets
            if char == delimiter_stack[-1]:  # find matching delimiter
                delimiter_stack.pop()
            elif char in delimiters:
                delimiter_stack.append(delimiters[char])  # nested brackets
            elif char == '\\':
                escaped = True  # so escaped quotes can be ignored
            else:
                kv_pairs[-1][-1].append(char)  # inside quotes, add verbatim
        elif char == '\\':
            escaped = True
        elif char in delimiters:
            delimiter_stack.append(delimiters[char])  # brackets or quotes
        elif (sep is None and char.isspace()) or char == sep:
            if kv_pairs == [[[]]]:  # empty, beginning of string
                continue
            elif kv_pairs[-1][-1] == []:
                continue
            else:
                kv_pairs.append([[]])
        elif char == '=':
            if kv_pairs[-1] == [[]]:
                del kv_pairs[-1]
            kv_pairs[-1].append([])  # value
        else:
            kv_pairs[-1][-1].append(char)

    kv_dict = {}

    for kv_pair in kv_pairs:
        if len(kv_pair) == 0:  # empty line
            continue
        elif len(kv_pair) == 1:  # default to True
            key, value = ''.join(kv_pair[0]), 'T'
        else:  # Smush anything else with kv-splitter '=' between them
            key, value = ''.join(kv_pair[0]), '='.join(
                ''.join(x) for x in kv_pair[1:])

        if key.lower() not in UNPROCESSED_KEYS:
            # Try to convert to (arrays of) floats, ints
            split_value = re.findall(r'[^\s,]+', value)
            try:
                try:
                    numvalue = np.array(split_value, dtype=int)
                except (ValueError, OverflowError):
                    # don't catch errors here so it falls through to bool
                    numvalue = np.array(split_value, dtype=float)
                if len(numvalue) == 1:
                    numvalue = numvalue[0]  # Only one number
                elif len(numvalue) == 9:
                    # special case: 3x3 matrix, fortran ordering
                    numvalue = np.array(numvalue).reshape((3, 3), order='F')
                value = numvalue
            except (ValueError, OverflowError):
                pass  # value is unchanged

            # Parse boolean values: 'T' -> True, 'F' -> False,
            #                       'T T F' -> [True, True, False]
            if isinstance(value, str):
                str_to_bool = {'T': True, 'F': False}

                try:
                    boolvalue = [str_to_bool[vpart] for vpart in
                                 re.findall(r'[^\s,]+', value)]
                    if len(boolvalue) == 1:
                        value = boolvalue[0]
                    else:
                        value = boolvalue
                except KeyError:
                    pass  # value is unchanged

        kv_dict[key] = value

    return kv_dict


if __name__ == '__main__':
    # Complex properties line. Keys and values that break with a regex parser.
    # see https://gitlab.com/ase/ase/issues/53 for more info

    complex_xyz_string = (
        ' '  # start with a separator
        'str=astring '
        'quot="quoted value" '
        u'quote_special="a_to_Z_$%%^&*\xfc\u2615" '
        r'escaped_quote="esc\"aped" '
        'true_value '
        'false_value = F '
        'integer=22 '
        'floating=1.1 '
        'int_array={1 2 3} '
        'float_array="3.3 4.4" '
        'a3x3_array="1 4 7 2 5 8 3 6 9" '  # fortran ordering
        'Lattice="  4.3  0.0 0.0 0.0  3.3 0.0 0.0 0.0  7.0 " '  # spaces in array
        'scientific_float=1.2e7 '
        'scientific_float_2=5e-6 '
        'scientific_float_array="1.2 2.2e3 4e1 3.3e-1 2e-2" '
        'not_array="1.2 3.4 text" '
        'nested_brackets=[[1,2],[3,4]] '  # gets flattented if not 3x3
        'bool_array={T F T F} '
        'bool_array_2=" T, F, T " '  # leading spaces
        'not_bool_array=[T F S] '
        # read and write
        u'\xfcnicode_key=val\xfce '
        u'unquoted_special_value=a_to_Z_$%%^&*\xfc\u2615 '
        '2body=33.3 '
        'hyphen-ated '
        # parse only
        'many_other_quotes=({[4 8 12]}) '
        'comma_separated="7, 4, -1" '
        'bool_array_commas=[T, T, F, T] '
        'Properties=species:S:1:pos:R:3 '
        'multiple_separators       '
        'double_equals=abc=xyz '
        'trailing'
    )

    expected_dict = {
        'str': 'astring',
        'quot': "quoted value",
        'quote_special': u"a_to_Z_$%%^&*\xfc\u2615",
        'escaped_quote': r"esc\"aped",
        'true_value': True,
        'false_value': False,
        'integer': 22,
        'floating': 1.1,
        'int_array': np.array([1, 2, 3]),
        'float_array': np.array([3.3, 4.4]),
        'a3x3_array': np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        'Lattice': np.array([[4.3, 0.0, 0.0], [0.0, 3.3, 0.0], [0.0, 0.0, 7.0]]),
        'scientific_float': 1.2e7,
        'scientific_float_2': 5e-6,
        'scientific_float_array': np.array([1.2, 2200, 40, 0.33, 0.02]),
        'not_array': "1.2 3.4 text",
        'nested_brackets': np.array([1, 2, 3, 4]),
        'bool_array': np.array([True, False, True, False]),
        'bool_array_2': np.array([True, False, True]),
        'not_bool_array': 'T F S',
        u'\xfcnicode_key': u'val\xfce',
        'unquoted_special_value': u'a_to_Z_$%%^&*\xfc\u2615',
        '2body': 33.3,
        'hyphen-ated': True,
        'many_other_quotes': np.array([4, 8, 12]),
        'comma_separated': np.array([7, 4, -1]),
        'bool_array_commas': np.array([True, True, False, True]),
        'Properties': 'species:S:1:pos:R:3',
        'multiple_separators': True,
        'double_equals': 'abc=xyz',
        'trailing': True
    }

    parsed_dict = key_val_str_to_dict(complex_xyz_string)
    np.testing.assert_equal(parsed_dict, expected_dict)
    print(parsed_dict)
