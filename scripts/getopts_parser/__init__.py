#!/usr/bin/env python3
# encoding=utf-8
"""!
@package getopts_parser
@brief Directly extract options from a command line into a dict. All option's values are cast according to a dictionary
which specifies options short and long names, their aliases and their type and their default value.

@section modules Modules and key functions
- @ref getopts_parser.py: Contain all code used to parse a command line.
    - getopts : The main function. Handles error.
    - getopts_parser : Same as getopts but without any error handling

- @ref tests.py: Provides a number of test function to verify the correct functioning the getopts_parser.py.
    - test_all : Call all test functions

@section usage_examples Usage Examples
Basic usage:
@code{.py}
import mypackage
mypackage.init_function()
# Output: This is an init function!

from mypackage import module1
module1.module1_function()
@endcode

Advanced usage:
@code{.py}
options_ = {
    'Alpha=': ("a:", (None, int)),
    'Beta=': ("b:", (None, list)),
    'Gamma': ("g", (True, False)),
    'Delta=': None,
    'Epsilon=': (None, ("Star", None)),
    'Zeta': ("z", (True, False)),
    'Eta': None,
    'Theta': None,
    'Iota': None,
    ('Kappa', 'Kap'): None,

    'Lambda': None,
    'Mu': None,
}

val = getopts_parser.getopts("1 --Eta --Iota -b Test2 --Kap" -g,
                             options_, fill_with_default_values=True)

#
# val["Alpha"] == 1
# val["Beta"] == list("Test2")
# val["Gamma"] is True
# val["Epsilon"] == "Star"
# val["Kappa"] is True
# val["Gamma"] is False
# val["Zeta"] is True
print("\tSuccess\n")
@endcode

@section dependencies Dependencies
- getopt : Parse a command line.
- pytest : Make test.py compatible with github workflow.
"""
from . import getopts_parser
from . import tests
