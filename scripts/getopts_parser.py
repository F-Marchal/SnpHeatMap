import getopt


class GetoptDigestionError(getopt.GetoptError):
    pass


class GetoptKeyError(getopt.GetoptError):
    pass

class GetoptParsingError(getopt.GetoptError):
    pass

def getopts_parser_simple_value(value_restriction: tuple[any, any], use_default=False) -> any:
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


def getopts_parser_complex_value(value_restriction, value, use_default=False) -> any:
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
        if callable(value_restriction[0]):
            return value_restriction[0]()
        else:
            return value_restriction[0]

    elif value_restriction[1] is None:
        return value
    else:
        return value_restriction[1](value)


def getopts_digester_check_item_endings(list_of_values: list[str], last_char: str = ":") -> bool or None:
    have_value = None

    for i, values in enumerate(list_of_values):
        if values == "":
            raise ValueError(f"Can not use empty strings as option")

        if have_value is None:
            have_value = values[-1] == last_char

        elif have_value is True and values[-1] == last_char:
            list_of_values[i] = values[:-1]

        elif have_value is False and values[-1] != last_char:
            pass

        else:
            raise ValueError(f"Conflict found. Some value(s) end by {last_char} but other(s) do not.")

    return have_value


def getopts_digest_available_options(getopts_options):
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

        else:
            # Erase duplicates
            short_keys = list(set(short_keys))

        short_string += "".join(short_keys)
        short_keys = ["-" + str(keys) for keys in short_keys]  # Add suffix

        # Check for incoherent options
        try:
            long_keys_ask_for_values = getopts_digester_check_item_endings(long_keys, last_char="=")
            short_keys_ask_for_values = getopts_digester_check_item_endings(short_keys, last_char=":")

            if short_keys_ask_for_values is not None and short_keys_ask_for_values != long_keys_ask_for_values:
                raise ValueError("Incoherent long keys / short keys. Short keys and long keys does not agree on"
                                 "whether they're waiting for a value")

        except ValueError as ValEr:
            raise GetoptDigestionError(msg=f"Incoherent options in the getopts' dict (key='{long_keys}', "
                                           f"values={short_keys_and_defaults}) :\n" + str(ValEr))

        # Digest default
        if defaults is None:
            if long_keys_ask_for_values:
                defaults = (None, None)
            
            else:
                defaults = (False, True)

        # Save results
        option_dict[main_key] = defaults

        # Store available keys
        for keys in [*long_keys, *short_keys]:
            if keys in complex_keys or keys in boolean_keys:
                raise GetoptKeyError(msg=f"Redundant option : {keys}", opt=keys)

            if long_keys_ask_for_values is True:
                complex_keys[keys] = main_key
            else:
                boolean_keys[keys] = main_key

    return option_dict, boolean_keys, complex_keys, short_string, long_list


def getopts_retrieve_options(argv, short_string, long_list, main_options):
    try:
        # Obtain values from argv
        opts, unparsed = getopt.getopt(argv, short_string, long_list)

        # When argument are passed
        if len(opts) == 0 and len(unparsed) >= 1:
            i = 0

            # Add an option to arguments
            #     |---------- Length verification ----------|  |- Option that accept value ? -||--- Not an option ---|
            while i < len(unparsed) and i < len(main_options) and main_options[i][-1] == "=" and unparsed[i][0] != "-":
                unparsed.insert(i * 2, "--" + main_options[i][:-1])
                i += 1
            print("la")
            opts, unparsed = getopt.getopt(unparsed, short_string, long_list)

        # Some options can not be parsed
        if unparsed:
            er_op = "".join(unparsed)
            raise GetoptParsingError(f"Unparsable option (Look for missing '-') : {er_op}", opt=er_op)

    except getopt.GetoptError as GetE:
        raise GetoptParsingError(msg=str(GetE))

    return opts


def getopts_parser(argv: list[str], getopts_options: dict[str or tuple, None or tuple[any, any]],
                   *mandatory: str, help_message: str = None, force_default=True) -> dict:
    """!
    @brief Use @ref getopt.getopt to generate a dictionary for a function. All items can have default values and
    returned values can be trans typed (or passed inside a function)

    @param argv : list[str] => List of argument (strings). Usually sys.argv[1:].
    @note If the first values are not known option (option in @p getopts_options), the function will consider that
    those values correspond to the nth first option that require an argument in @p getopts_options.
    In below example, you can give two arguments. The first will be matched with "Alpha" and the last with "Beta".
    @code
    argv = ["5", "6"]
    getopts_options = {
        "Alpha=":               ("a:", None),
        "Beta=:                 ('b', (0, int)),
        "Gamma":                ("g", ("data/", str)),
        "Delta=":               ("d:", (0, lambda integer : int(integer) - 1)),
    }
    @endcode

    @param getopts_options : dict[str,None or tuple[any, any]] => dictionary that contain.
    This dict should have the following format : {complete_name: (shorten_name, (default_value, value_restriction) ) }
        - complete_name : Option called using --Complete_name. It's this name that is used to fill the returned dict.
        - shorten_name : Option called using -s.
        @note If an option require an argument, complete_name should end by an '=' AND shorten_name should end by a ':'.
        - (default_value, value_restriction):
            - Option that does not require an argument :
                - When None : None is replaced by (True, False)
                - default_value : The value used when this option is not used.
                - value_restriction : The value used when this option is used.

            - Option that require an argument :
                - None : Only for options in @p mandatory
                - default_value : The value used when this option is not used. (can be a callable)
                - value_restriction : None or a callable. Use it to turn your value in other type e.g. : int
    Example :
    @code
    getopts_options = {
        "Alpha=":               ("a:", None),
        "Beta":                 (None, None),
        "Gamma":                ("g", ("data/", str)),
        "Delta=":               ("d:", (0, lambda integer : int(integer) - 1)),
    }
    @endcode

    @param *mandatory : str => A list of option / argument that should be used at each time (e.g. file path in cut
    command)
    @param help_message : str = None => A message that will be displayed if the arguments are incorrect or if
    --help is used (any short key can be attributed to --help)
    @param force_default = True => All unspecified values are filled with theirs default value (see @p getopts_options)

    @note if a parameter is named "help" and that @p help_message is True,
    @return dict => list of parameter xtract from @p argv

    """

    option_dict, boolean_keys, complex_keys, short_string, long_list = getopts_digest_available_options(getopts_options)
    main_options = sorted(option_dict)

    options = getopts_retrieve_options(argv, short_string, long_list, main_options)

    for opt in options:
        if opt in boolean_keys:
            pass
    '''option_dict = {}            # Will store options

    for option, value in opts:

        if option in simple_keys:
            # python3 main.py -rf --> [("-r", "f")].
            # That not what we want. We want  :
            # python3 main.py -rf -->[("-r", ""), ("-f", "")]:
            # In order to achieve that "value" at the end of opts

            if value:
                opts.extend([("-" + other_options, "") for other_options in value])

            # Apply filter
            value_restriction, option_name = simple_keys[option]
            option_dict[option_name] = getopts_parser_simple_value(value_restriction)

        elif option in complexe_keys:
            value_restriction, option_name = complexe_keys[option]
            option_dict[option_name] = getopts_parser_complex_value(value_restriction, value)

        else:
            raise KeyError("Unknown option, can not proceed : " + option)

        if option_name in mandatory:
            mandatory.remove(option_name)  # option_name comes from precedent if / else

    if mandatory:
        raise KeyError(f"{len(mandatory)} argument missing : {', '.join(mandatory)}")

    if help_message and "help" in option_dict:
        del option_dict["help"]

    return option_dict'''


__getops__  = {
    "Ug=": None,
    "Opt=": None,
    ("LongKey1", "LongKey2"): (None, None),
}

# Multiple / mono long keys
# Multiple / mono small keys
# No small key
import sys

getopts_parser(sys.argv[1:], __getops__)
