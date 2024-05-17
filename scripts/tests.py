#!/usr/bin/env python3
# encoding=utf-8
"""! @brief Script to use in order to test "snp_analyser" functionalities.
 @file tests.py
 @section libs Libraries / Modules
  - snp_analyser
  - pytest
 @section authors Author(s)
  - Created by Marchal Florent on 17/5/2024 .
"""
try:
    import scripts.snp_analyser as snp
    import pytest
except ModuleNotFoundError as E:
    print(f"Module not found : {E}\n")

    if E.name == "snp_analyser":
        print(f"Make sure that 'snp_analyser' is in the same folder as {__file__}")

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

# def test_help_usage(): # Nothing to do


def test_greater_than_0_int_filter():
    """@brief Test greater_than_0_int_filter"""
    assert snp.filter_integer_greater_or_equal_to_0_ignore_0(0) is False
    assert snp.filter_integer_greater_or_equal_to_0_ignore_0(+2000) is True
    assert snp.filter_integer_greater_or_equal_to_0_ignore_0(1_000) is True
    assert snp.filter_integer_greater_or_equal_to_0_ignore_0("+2000") is True
    assert snp.filter_integer_greater_or_equal_to_0_ignore_0("20_00") is True

    with pytest.raises(ValueError):
        assert snp.filter_integer_greater_or_equal_to_0_ignore_0(-1000)
        assert snp.filter_integer_greater_or_equal_to_0_ignore_0("twenty")
        assert snp.filter_integer_greater_or_equal_to_0_ignore_0("0.1")
        assert snp.filter_integer_greater_or_equal_to_0_ignore_0("1.0")


def test_compile_gene_snp():
    """@brief Test compile_gene_snp"""
    genes1 = {
        'gene0': 0,
        'gene1': 3,
        'gene2': 4,
        'gene3': 4,
        'gene4': 1,
        'gene5': 4,
        'gene6': 2,
        'gene7': 3,
        'gene8': 2,
        'gene9': 1
    }

    genes1_a = {
        'gene0': 0,
        'gene1': 3,
        'gene2': 4,
        'gene3': 4,
        'gene4': 1,
    }
    genes1_b = {
        'gene5': 4,
        'gene6': 2,
        'gene7': 3,
        'gene8': 2,
        'gene9': 1
    }

    genes1_b_bis = {
        'gene0': 0,
        'gene1': 2,
        'gene2': 1,
        'gene3': 4,
        'gene4': 3,
    }

    genes_c_bis = {
        'gene0': 0,
        'gene1': 2,
        'gene2': 1,
        'gene3': 4,
        'gene4': 3,
    }

    none_result = {0: {'None': 1}, 3: {'None': 2}, 4: {'None': 3}, 1: {'None': 2}, 2: {'None': 2}}
    john_smith_result = {0: {'John Smith': 1}, 3: {'John Smith': 2}, 4: {'John Smith': 3}, 1: {'John Smith': 2}, 2:
                         {'John Smith': 2}}
    a_b_result = {0: {'a': 1}, 3: {'a': 1, 'b': 1}, 4: {'a': 2, 'b': 1}, 1: {'a': 1, 'b': 1}, 2: {'b': 2}}
    a_b_result2 = {0: {'a': 1, 'b': 1}, 3: {'a': 1, 'b': 2}, 4: {'a': 2, 'b': 2}, 1: {'a': 1, 'b': 2}, 2: {'b': 3}}
    a_c_b_result = {0: {'a': 1, 'b': 1, 'c': 1}, 3: {'a': 1, 'b': 2, 'c': 1}, 4: {'a': 2, 'b': 2, 'c': 1},
                    1: {'a': 1, 'b': 2, 'c': 1}, 2: {'b': 3, 'c': 1}}

    assert snp.compile_gene_snp(genes1) == none_result
    assert snp.compile_gene_snp(genes1, group="John Smith") == john_smith_result

    temp = snp.compile_gene_snp(genes1_a, group="a")
    temp = snp.compile_gene_snp(genes1_b, group="b", dict_of_number=temp)
    assert temp == a_b_result
    assert snp.compile_gene_snp(genes1_b_bis, group="b", dict_of_number=temp) == a_b_result2
    assert snp.compile_gene_snp(genes_c_bis, group="c", dict_of_number=temp) == a_c_b_result


def test_make_data_matrix():
    """@brief Test make_data_matrix"""
    a_c_b_result = {3: {'a': 3, 'b': 2, 'c': 1}, 4: {'a': 8, 'b': 2, 'c': 1},
                    1: {'a': 7, 'b': 2, 'c': 1}, 2: {'b': 8, 'c': 6}, 6: {'b': 3}, 12: {'b': 6}, 8: {'c': 7}}
    # assert snp.make_data_matrix(a_c_b_result, "a") == ([[7, 0, 3, 8]], [1, 2, 3, 4])
    print()
    assert snp.make_data_matrix(a_c_b_result, "b", simplified=False) == ([[2, 8, 2, 2, 0, 3, 0, 0, 0, 0, 0, 6]],
                                                                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    assert snp.make_data_matrix(a_c_b_result, "b", simplified=True) == ([[2, 8, 2, 2, 3, 0, 6]],
                                                                        [1, 2, 3, 4, 6, 8, 12])
    assert snp.make_data_matrix(a_c_b_result, "c", simplified=False) == ([[1, 6, 1, 1, 0, 0, 0, 7, 0, 0, 0, 0]],
                                                                          [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    assert snp.make_data_matrix(a_c_b_result, "c", simplified=True) == ([[1, 6, 1, 1, 0, 7, 0]],
                                                                        [1, 2, 3, 4, 6, 8, 12])

    assert snp.make_data_matrix(a_c_b_result, "a", simplified=False) == ([[7, 0, 3, 8, 0, 0, 0, 0, 0, 0, 0, 0]],
                                                                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    assert snp.make_data_matrix(a_c_b_result, "a", simplified=True) == ([[7, 0, 3, 8, 0, 0, 0]],
                                                                        [1, 2, 3, 4, 6, 8, 12])

    assert snp.make_data_matrix(a_c_b_result, "b", 'c', simplified=False) == ([
                                                                                 [2, 8, 2, 2, 0, 3, 0, 0, 0, 0, 0, 6],
                                                                                 [1, 6, 1, 1, 0, 0, 0, 7, 0, 0, 0, 0]
                                                                             ],
                                                                             [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    assert snp.make_data_matrix(a_c_b_result, "a", 'c', simplified=True) == ([[7, 0, 3, 8, 0, 0, 0],
                                                                             [1, 6, 1, 1, 0, 7, 0]],
                                                                            [1, 2, 3, 4, 6, 8, 12])
    assert snp.make_data_matrix(a_c_b_result, "b", 'c', "a", simplified=False) == ([
                                                                                   [2, 8, 2, 2, 0, 3, 0, 0, 0, 0, 0, 6],
                                                                                   [1, 6, 1, 1, 0, 0, 0, 7, 0, 0, 0, 0],
                                                                                   [7, 0, 3, 8, 0, 0, 0, 0, 0, 0, 0, 0]
                                                                                  ],
                                                                                  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    assert snp.make_data_matrix(a_c_b_result, "a", 'c', "b", simplified=True) == ([[7, 0, 3, 8, 0, 0, 0],
                                                                                   [1, 6, 1, 1, 0, 7, 0],
                                                                                   [2, 8, 2, 2, 3, 0, 6]
                                                                                   ],
                                                                                  [1, 2, 3, 4, 6, 8, 12])
    assert snp.make_data_matrix(a_c_b_result, "a", simplified=False, max_length=5) == ([[7, 0, 3, 8, 0, 0]],
                                                                                       [1, 2, 3, 4, 5, 6])
    assert snp.make_data_matrix(a_c_b_result, "a", simplified=True, max_length=5) == ([[7, 0, 3, 8, 0]],
                                                                                      [1, 2, 3, 4, 6])


def test_generate_cumulative_list():
    "list_of_numbers: list[int] or list[float], reversed_=False"
    assert snp.generate_cumulative_list([3, 35, 2, 1, 3, 4]) == [3, 38, 40, 41, 44, 48]
    assert snp.generate_cumulative_list([3, 35, 2, 1, 3, 4], reversed_=True) == [48, 45, 10, 8, 7, 4]
    assert snp.generate_cumulative_list([0, 0, 2, 1, 0.35, 0]) == [0, 0, 2, 3, 3.35, 3.35]
    assert snp.generate_cumulative_list([0, 0, 2, 1, 0.35, 0], reversed_=True) == [3.35, 3.35, 3.35, 1.35, 0.35, 0]


def test_all():
    """
    @brief Call all test functions
    """
    test_compile_gene_snp()
    test_greater_than_0_int_filter()
    test_make_data_matrix()
    test_generate_cumulative_list()


if __name__ == '__main__':
    test_all()

