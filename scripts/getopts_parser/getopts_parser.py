#!/usr/bin/env python3
# encoding=utf-8
"""!Directly extract options from a command line into a dict. All option's values are cast according to a dictionary
which specifies options short and long names, their aliases and their type and their default value.

@see getopts
@note You can read test.py if you need some usage example of getopts

@file getopts_parser.py
@section libs Librairies/Modules
- getopt
@section authors Author(s)
- Created by Marchal Florent on 16/05/2024 .
"""
import getopt


class GetoptsDigestionError(getopt.GetoptError):
    """!
    @brief Error raised when an occur during the digestion of the dictionary that contain options accepted by a getopts
    """


class GetoptsOptionError(getopt.GetoptError):
    """!
    @brief Error raised when there is an error with an option related to a getopts
    """


class GetoptsParsingError(getopt.GetoptError):
    """!
    @brief Error raised when there is an error during the parsing of a command-line.
    """


def getopts_parser_boolean_option(value_restriction: tuple[any, any] or None, use_default=False) -> any:
    """!
    @brief Internal function used by  @ref getopts_parser.

    @param value_restriction : tuple[any, any] => A tuple of two value. (Default value, normal value).
    (False, True) when None
    @param use_default = False => Do the default value is returned

     @return any =>
        - True if value_restriction is None
        - Default value if @p value_restriction is True (@p value_restriction [1])
        - Normal value if @p value_restriction is False (@p value_restriction [0])
    """
    if value_restriction is None:
        value_restriction = (False, True)

    if use_default:
        return value_restriction[0]

    else:
        return value_restriction[1]


def getopts_parser_complex_option(value_restriction, value, use_default=False) -> any:
    """!
    @brief Internal function used by  @ref getopts_parser. Use it to apply a number of restriction on a value.

    Parameters :
        @param value_restriction : tuple[any, callable] => > A tuple of two value. (Default value, normal value).
         (False, True) when None
         - if value_restriction[0] is callable : default = result of the function (no argument is given)
         - if value_restriction[0] is not callable : default = value_restriction[0]
         - if value_restriction[1] callable : is applied to @p value, the result is used as return value. If None
         @p value is returned with no change
        @param value => Value that can be returned. This value pass inside @p value_restriction [1]
        @param use_default = False => Do the default value is returned

     @return any => @p value transformed by  @p value_restriction [1]
    """
    if value_restriction is None:
        value_restriction = (None, None)

    if use_default:
        # Default value is callable
        if callable(value_restriction[0]):
            return value_restriction[0]()

        # Default value is not callable
        else:
            return value_restriction[0]

    # Cast @p value using value_restriction[1]
    elif value_restriction[1] is None:
        return value
    else:
        return value_restriction[1](value)


def getopts_digester_check_item_endings(list_of_values: list[str], last_char: str = ":",) -> bool or None:
    """!
    @brief Check if the last char of each string in @p list_of_values is equal to last_char.
    if some od themes are equals and others are not, an error is raised.

    @param list_of_values : list[str] => A list of string
    @param last_char : str = ":" => A character

    @return bool or None =>
    - True : All items in @p list_of_values end by @p last_char
    - False : No item in @p list_of_values end by @p last_char
    - None : @p list_of_values is empty
    - ValueError :  Some value(s) end by @p last_char but other(s) do(es) not."
    """

    have_value = None

    for i, values in enumerate(list_of_values):
        # Do the first item end by @p last_char ?
        if have_value is None:
            have_value = values[-1] == last_char

        # Test if this item end by the same symbol as the first one
        elif have_value is True and values[-1] == last_char:
            pass

        elif have_value is False and values[-1] != last_char:
            pass

        else:
            raise ValueError(f"Conflict found. Some value(s) end by '{last_char}' but other(s) do(es) not.")

    return have_value


def getopts_digest_available_options(getopts_options: dict[str, tuple[any, any]],
                                     dict_of_default_value: dict[str, any] = None) \
        -> tuple[dict[str, tuple[any, any]], dict[str, str], dict[str, str], str, list[str]]:
    """!
    @brief Transform @p getopts_options into usable instruction for @ref getopts_retrieve_options and assure that
    option's aliases are coherent.


     @param getopts_options : dict[str,tuple[any,any]] =>  A dictionary that follows one of the following structure :
    @code
    {name: None}
    {tuple_of_name: None}
    {name: (short_name, None)}
    {tuple_of_name: (tuple_of_short_name, None)}

    {name: (tuple_of_short_name, default_value)}
    {name: (short_name, default_value)}
    {tuple_of_name: (tuple_of_short_name, default_value)}
    {tuple_of_name: (short_name, default_value)}

    {name: (tuple_of_short_name, (default_value, cast))}
    {name: (short_name, (default_value, cast))}
    {tuple_of_name: (tuple_of_short_name, (default_value, cast))}
    {tuple_of_name: (short_name, (default_value, cast))}
    @endcode
    @param dict_of_default_value : dict[str,any] = None => Facultative, A dictionary that will be filled using options's
    default values.


    @return tuple[dict[str, tuple[any, any]], dict[str, str], dict[str, str], str, list[str]] =>
    - option_dict : dict[str, tuple[any, any] => Contain each option "main_name" associated with theirs
        default and cast options
    - boolean_keys : dict[str, str] => Contain all options aliases related to option that can handle only two state
    - complex_keys : dict[str, str] => Contain all options aliases related to option that can handle more than two state
    - short_string : str => string that correspond to shortopts in  @ref getopt.getopt
    - long_list : list[str] => list of string that correspond to longopts in @ref getopt.getopt

    """

    boolean_keys = {}
    complex_keys = {}
    short_string = ""
    long_list = []
    option_dict = {}

    for long_keys, short_keys_and_defaults in getopts_options.items():
        # Digest long keys
        if isinstance(long_keys, tuple):
            if not long_keys:
                continue
            # Erase duplicates
            main_key = long_keys[0]
            long_keys = list(set(long_keys))
        
        else:
            main_key = long_keys
            long_keys = [long_keys]

        long_list.extend(long_keys)
        long_keys = ["--" + str(keys) for keys in long_keys]   # Add suffix

        # Split short keys and defaults
        if short_keys_and_defaults is None:
            short_keys = None
            defaults = None

        else:
            # Unpack short_keys_and_defaults
            short_keys, defaults = short_keys_and_defaults

        # Digest short keys
        if short_keys is None:
            short_keys = []

        elif isinstance(short_keys, (tuple, list, set)):
            # Erase duplicates
            short_keys = list(set(short_keys))

        else:
            short_keys = [short_keys]

        short_string += "".join(short_keys)

        # add suffix
        temp = []  # Add suffix
        for keys in short_keys:
            if not keys:
                continue
            if (len(keys) == 1 and ":" in keys) or (len(keys) == 2 and keys[1] != ":") or "-" in keys or "=" in keys:
                raise GetoptsDigestionError("Short keys should contain only one car and can not contain"
                                            f"'=' or '-'. Got {keys}", opt=keys)
            temp.append("-" + keys)

        short_keys = temp

        del temp

        # Check for incoherent options
        try:
            long_keys_ask_for_values = getopts_digester_check_item_endings(long_keys, last_char="=")
            short_keys_ask_for_values = getopts_digester_check_item_endings(short_keys, last_char=":")

            if short_keys_ask_for_values is not None and short_keys_ask_for_values != long_keys_ask_for_values:
                raise ValueError("Incoherent long keys / short keys. Short keys and long keys does not agree on "
                                 "whether they're waiting for a value")

        except ValueError as ValEr:
            raise GetoptsDigestionError(msg=f"Incoherent options in the getopts' dict (key='{long_keys}', "
                                            f"values={short_keys}) :\n" + str(ValEr))

        # Digest default
        if defaults is None or defaults == tuple():
            if long_keys_ask_for_values:
                defaults = (None, None)
            
            else:
                defaults = (False, True)

        elif not isinstance(defaults, tuple):
            defaults = (None, defaults)

        # Save results
        option_dict[main_key] = defaults

        # Store available keys
        for keys in [*long_keys, *short_keys]:
            if keys in complex_keys or keys in boolean_keys:
                raise GetoptsOptionError(msg=f"Redundant option : {keys}", opt=keys)

            if long_keys_ask_for_values is True:
                complex_keys[keys] = main_key

                if dict_of_default_value is not None:
                    dict_of_default_value[main_key[:-1]] = getopts_parser_complex_option(defaults, None, True)

            else:
                boolean_keys[keys] = main_key

                if dict_of_default_value is not None:
                    dict_of_default_value[main_key] = getopts_parser_boolean_option(defaults, True)

    return option_dict, boolean_keys, complex_keys, short_string, long_list


def getopts_retrieve_options(argv: list[str], short_string: str, long_list: list[str], main_options: list[str]) \
        -> list[tuple[str, str]]:
    """!
    @brief Use a list of string from a command-line to extract options and theirs values.

    @param argv : list[str] => A list of string (Same format as sys.argv[1:])
    @param short_string : str => string that correspond to shortopts in  @ref getopt.getopt
    @param long_list : list[str] => list of string that correspond to longopts in @ref getopt.getopt
    @param main_options : list[str] => An ordered list of all options. Is used to determine the
    correspondence between arguments and parameters

    @return list[Tuple[str, str]] => A list composed of tuple that contain option and their values.

    """
    try:
        # Obtain values from argv
        opts, unparsed = getopt.getopt(argv, short_string, long_list)

        # When argument are passed
        if len(opts) == 0 and len(unparsed) >= 1:
            i = 0

            # Add an option to arguments
            #     |---------- Length verification ----------| |- Option that accept value ? -|
            while (i < len(unparsed) and i < len(main_options) and main_options[i][-1] == "="
                   and unparsed[i * 2][0] != "-"):
                #       |--- Not an option ---|

                unparsed.insert(i * 2, "--" + main_options[i][:-1])
                i += 1

            opts, unparsed = getopt.getopt(unparsed, short_string, long_list)

        # Some options can not be parsed
        if unparsed:
            er_op = " ".join(unparsed)
            if len(opts) == 0:
                raise GetoptsParsingError(f"Unparsable option (Look for missing '-') : {er_op}", opt=er_op)

            else:
                func_op = " ".join([f"{key} {value}" if value else f"{key}" for key, value in opts])
                raise GetoptsParsingError(f"Unable to parse all option. (Look for missing '-' and parameters"
                                          f"which have a parameter entered even if it should not)"
                                          f"\nDid : {func_op}"
                                          f"\nMissing : {er_op}", opt=er_op)
    except getopt.GetoptError as GetE:
        raise GetoptsParsingError(msg=str(GetE))

    return opts


def getopts_parser(argv: list[str] or str, getopts_options: dict[str, tuple[any, any]],
                   *mandatory: str, fill_with_default_values=True) -> dict:
    """!
    @brief Internal version of @ref getopts. You can use this function if you don't need @ref getopts ' overcoat.
    (Display help message and error handling)

    @see getopts


    @param argv : list[str] or str => List of argument (strings). Usually sys.argv[1:] or str command line.
    @param getopts_options : dict[str, tuple[any, any]] => A dictionary that contain options.
    @param *mandatory : str =>a list of options whose value must be entered
    @param fill_with_default_values = True => Every unused options are added to the final dict using defaults values

    @return dict => A dictionary that contain options.

    """
    final_values = {}
    mandatory = set(mandatory)

    # str compatibility
    if isinstance(argv, str):
        argv = argv.split(" ")

    # digest @p getopts_options
    if fill_with_default_values:
        option_dict, boolean_keys, complex_keys, short_string, long_list \
            = getopts_digest_available_options(getopts_options, final_values)
    else:
        option_dict, boolean_keys, complex_keys, short_string, long_list \
            = getopts_digest_available_options(getopts_options)

    # extract options from argv
    options = getopts_retrieve_options(argv, short_string, long_list, main_options=list(option_dict))

    # Fill final_values
    for opt in options:
        opt, value = opt

        if opt in boolean_keys:
            # "-opt"
            if value:
                # Sometime, when there is multiple chars after a '-' (e.g. -atcg) some letters are considered as value :
                # --> -atcg <=W ('-a', 'tcg')
                options.extend([("-" + char, "") for char in value])

            # Apply value restriction
            parent = boolean_keys[opt]
            value = getopts_parser_boolean_option(option_dict[parent])

            if parent in mandatory:
                mandatory.remove(parent)
            
        else:
            if value and value[0] == "-" and value in complex_keys or value in boolean_keys:
                # Avoid that an option is taken as value by another one
                # Still allow negative numbers
                raise GetoptsOptionError(f"The option '{opt}' is followed by a parameter '{value}' instead of "
                                         f"a value.")

            # Apply value restriction
            if opt[0:2] == "--":
                parent = complex_keys[opt + "="]

            else:
                parent = complex_keys[opt + ":"]

            value = getopts_parser_complex_option(option_dict[parent], value)
            parent = parent[:-1]

            if parent in mandatory:
                mandatory.remove(parent)

        # Save result
        final_values[parent] = value

    # Verify that all mandatory values are in the dict
    if mandatory:
        opt = " ".join(mandatory)
        raise GetoptsParsingError(f"At least one option is missing : {opt}", opt=opt)

    return final_values


def getopts(argv: list[str] or str, getopts_options: dict[str or tuple, None or tuple[any, any] or any],
            *mandatory: str, help_message: str = None, fill_with_default_values=True,
            raise_errors=False, help_options: str or tuple[str] = ("help", )) -> dict[str, any] or int:
    """!
    @brief Directly extract options from a command line into a dict. All option's values are cast according to a
    dictionary @p getopts_options which specifies options short and long names, their aliases and their type and
    their default value.

    This function use @ref getopt.getopt.

    @param argv : list[str] or str => [description]
    @param getopts_options : dict[str, tuple[any, any]] =>
    - Structure for this dict :
    @code
    {name: None}
    {tuple_of_name: None}
    {name: (short_name, None)}
    {tuple_of_name: (tuple_of_short_name, None)}

    {name: (tuple_of_short_name, default_value)}
    {name: (short_name, default_value)}
    {tuple_of_name: (tuple_of_short_name, default_value)}
    {tuple_of_name: (short_name, default_value)}

    {name: (tuple_of_short_name, (default_value, cast))}
    {name: (short_name, (default_value, cast))}
    {tuple_of_name: (tuple_of_short_name, (default_value, cast))}
    {tuple_of_name: (short_name, (default_value, cast))}

    @endcode
    - Option that accept values should have a "=" at the end of each name and a ":" at the end of each short_name.
    - Option that require a value : (e.g. --file path) :
        - default_value can be of any type, if callable the result of the function is used
        - The cast can be any function. if None, nothing happen

    - Option that does not require any value : (e.g. --help) :
        - The default value is False
        - The "on" value is True
        - If you specie a "cast" it will be used as the "on" value

    @note When an option has multiple names, the first one in the tuple is selected as the "main_name". The "main_name"
    is the key used inside the returned dict.
    @note short_names can be composed by only one character (plus ":" if needed)
    @warning do not use "-" in your options names

    @param *mandatory : str => a list of options whose value must be entered
    @param help_message : str = None => A message displayed when an error is encountered or when "help" is triggered.
    @param fill_with_default_values = True => When True, the returned dict contain all options
    @param raise_errors = False => Do caught errors are raised ?
    @param help_options : str or tuple[str] = ("help", ) => A tuple of option name that trigger @p help_message.
    Those option are always removed from the returned dictionary.

    @return dict[str, any] or int =>
        -  dict[str, any] A dictionary that contain values related to options triggered by the command-line.
        - Only when an error occur and @p raise_errors is False or when help message is displayed.
            - 1 = An errors been caught
            - 2 = help message was required.
    """

    try:
        vals = getopts_parser(argv, getopts_options, *mandatory,
                              fill_with_default_values=fill_with_default_values)

    except getopt.GetoptError as E:
        print("An error occured during the parsing of your command-line arguments :"
              f"\n{E.msg}")

        if help_message:
            print(f"\nHere some help regarding the use of this script : \n{help_message}")

        if raise_errors is True:
            raise E

        return 1

    if not isinstance(help_options, tuple):
        help_options = [help_options]

    for help_opt in help_options:
        if help_opt not in vals:
            continue

        if vals[help_opt]:
            print(help_message)
            return 2

        else:
            del vals[help_opt]

    return vals
