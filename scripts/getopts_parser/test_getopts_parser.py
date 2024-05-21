#!/usr/bin/env python3
# encoding=utf-8
"""! @file test_getopts_parser.py
 @brief Script to use in order to test "getopts_parser" functionalities.
 @file test_snp_analyser.py
 @section libs Libraries / Modules
  - getopts_parser
  - pytest
 @section authors Author(s)
  - Created by Marchal Florent on 16/5/2024 .
"""
__author__ = "Marchal Florent"
__credits__ = ["Marchal Florent"]

try:
    from . import getopts_parser
    import pytest

except ModuleNotFoundError as E:
    print(f"Module not found : {E}\n")

    if E.name == "getopts_parser":
        print(f"Make sure that 'getopts_parser' is in the same folder as {__file__}")

    else:
        print(
          f"Open a terminal and try : "
          f"\n\tpip install pytest"
          
          f"\n\nIf pip is not found, you can install it using : "
          f"\nOn linux or MacOs: "
          f"\n\tpython -m ensurepip --upgrade"
          f"\nOn Windows : "
          f"\n\tpy -m ensurepip --upgrade"
          )
    exit(1)


def test_simple_option():
    """!
    @brief Test simple functionality of getopts_parser
    """
    options_ = {
        'Alpha=': None,
        'Beta=': None,
        'Gamma': None,
        'Delta=': None,
        'Epsilon=': None,
        'Zeta': None,
        'Eta': None,
        'Theta': None,
        'Iota': None,
        'Kappa': None,
        'Lambda': None,
        'Mu': None,
    }

    print("Arguments and parameters")
    val = getopts_parser.getopts("Test1 Test2 --Eta --Iota", options_, fill_with_default_values=False)
    assert not isinstance(val, int)
    assert val["Alpha"] == "Test1"
    assert val["Beta"] == "Test2"
    assert val["Eta"] is True
    assert val["Iota"] is True
    print("\tSuccess\n")

    print("Arguments using parameters")
    val = getopts_parser.getopts("Test1 --Eta --Iota --Beta Test2", options_, fill_with_default_values=False)
    assert not isinstance(val, int)
    assert val["Alpha"] == "Test1"
    assert val["Beta"] == "Test2"
    print("\tSuccess\n")

    print("Parameters erasing arguments")
    val = getopts_parser.getopts("Test1 --Eta --Iota --Alpha Test2", options_, fill_with_default_values=False)
    assert not isinstance(val, int)
    assert val["Alpha"] == "Test2"
    assert "Beta" not in val
    print("\tSuccess\n")

    print("Defaults Values")
    val = getopts_parser.getopts("Test1 --Eta --Iota --Alpha Test2", options_, fill_with_default_values=True)
    assert not isinstance(val, int)
    assert val["Alpha"] == "Test2"
    assert val["Iota"] is True
    assert val["Lambda"] is False
    print("\tSuccess\n")


def test_aliases():
    """!
    @brief Test aliases functionality
    """
    options_ = {
        'Alpha=': ("a:", None),
        'Beta=': ("b:", None),
        'Gamma': ("g", None),
        'Delta=': None,
        'Epsilon=': None,
        'Zeta': None,
        'Eta': None,
        'Theta': None,
        'Iota': None,
        'Kappa': None,
        'Lambda': None,
        'Mu': None,
    }

    print("Mono aliases (long + short)")
    val = getopts_parser.getopts("-a Test1 --Beta Test2 -g", options_, fill_with_default_values=False)
    assert not isinstance(val, int)
    assert val["Alpha"] == "Test1"
    assert val["Beta"] == "Test2"
    assert val["Gamma"] is True
    print("\tSuccess\n")

    options_ = {
        'Alpha=': ("a:", None),
        'Beta=': ("", None), # No short key
        'Gamma': None,
        ('Delta=', "delta=", "Del="): None,
        'Epsilon=': (("E:", "e:", "w:"), None),
        'Zeta': None,
        ('Eta', "eta", "ta", "Eta"): (("t", "r", "r"), None),
        'Theta': None,
        'Iota': (("I", "i", "v"), None),
        'Kappa': None,
        'Lambda': None,
        'Mu': None,
    }

    print("Multiple aliases (long + short)")
    val = getopts_parser.getopts("-a Test1 --Beta Test2 -i --ta -e Test4 --Del Test5",
                                 options_, fill_with_default_values=False)
    assert not isinstance(val, int)
    assert val["Alpha"] == "Test1"
    assert val["Beta"] == "Test2"
    assert val["Iota"] is True
    assert val["Eta"] is True
    assert val["Epsilon"] == "Test4"
    assert val["Delta"] == "Test5"
    print("\tSuccess\n")

    print("Some error test :")
    print("Incoherence:")
    options_["Bozo"] = ("b=", None)
    with pytest.raises(getopts_parser.GetoptsDigestionError):
        getopts_parser.getopts("-a Test1 --Beta Test -b temp", options_, raise_errors=True)

    print("\nIncoherence 2 :")
    options_["Bozo"] = (("b:", "c"), None)
    with pytest.raises(getopts_parser.GetoptsDigestionError):
        getopts_parser.getopts("-a Test1 --Beta Test -b temp", options_, raise_errors=True)

    print("\nToo long short :")
    options_["Bozo"] = ("bb", None)
    with pytest.raises(getopts_parser.GetoptsDigestionError):
        getopts_parser.getopts("-a Test1 --Beta Test -b temp", options_, raise_errors=True)

    print("\nIllegal short :")
    options_["Bozo"] = ("-", None)
    with pytest.raises(getopts_parser.GetoptsDigestionError):
        getopts_parser.getopts("-a Test1 --Beta Test -b temp", options_, raise_errors=True)

    print("\nOverlaps :")
    options_["Bozo"] = (("b:", "c", "i"), None)
    with pytest.raises(getopts_parser.GetoptsDigestionError):
        getopts_parser.getopts("-a Test1 --Beta Test -b temp", options_, raise_errors=True)

    print("\nIncoherence 3 :")
    del options_["Bozo"]
    options_[("Bozo", "Bozo2=", "Bozo3")] = (("b", "c"), None)
    with pytest.raises(getopts_parser.GetoptsDigestionError):
        getopts_parser.getopts("-a Test1 --Beta Test -b temp", options_, raise_errors=True)

    print("\tSuccess\n")


def test_defaults_values():
    """!
    @brief Test default value and cast
    """
    options_ = {
        'Alpha=': ("a:", None),
        'Beta=': ("b:", None),
        'Gamma': ("g", (True, False)),
        'Delta=': None,
        'Epsilon=': (None, ("Star", None)),
        'Zeta': None,
        'Eta': None,
        'Theta': None,
        'Iota': None,
        'Kappa': None,
        'Lambda': None,
        'Mu': None,
    }

    print("Simple Default :")
    val = getopts_parser.getopts("Test1 --Eta --Iota --Beta Test2", options_, fill_with_default_values=True)

    assert not isinstance(val, int)
    assert val["Alpha"] == "Test1"
    assert val["Beta"] == "Test2"
    assert val["Gamma"] is True
    assert val["Epsilon"] == "Star"
    print("\tSuccess\n")

    options_ = {
        'Alpha=': ("a:", (None, int)),
        'Beta=': ("b:", (None, list)),
        'Gamma': ("g", (True, False)),
        'Delta=': None,
        'Epsilon=': (None, ("Star", None)),
        'Zeta': None,
        'Eta': None,
        'Theta': None,
        'Iota': None,
        'Kappa': None,
        'Lambda': None,
        'Mu': None,
    }

    print("Default and cast :")
    val = getopts_parser.getopts("1 --Eta --Iota --Beta Test2", options_, fill_with_default_values=True)

    assert not isinstance(val, int)
    assert val["Alpha"] == 1
    assert val["Beta"] == list("Test2")
    assert val["Gamma"] is True
    assert val["Epsilon"] == "Star"
    print("\tSuccess\n")

    print("Testing parameter without option")
    with pytest.raises(getopts_parser.GetoptsOptionError):
        getopts_parser.getopts("1 --Eta --Iota --Beta -g", options_, raise_errors=True)

    with pytest.raises(getopts_parser.GetoptsParsingError):
        getopts_parser.getopts("1 --Eta --Iota --Beta", options_, raise_errors=True)
    print("\tSuccess\n")


def test_all():
    """!
    @brief Call each test functions
    """
    test_simple_option()
    test_aliases()
    test_defaults_values()


if __name__ == "__main__":
    test_all()
