"""!
@file snp_analyser.py
@brief
Create a number of chart related to snp (simple nucleotide polymorphism) analysis

@subsection How to use this file
Python3 [Column that contain gene's names] [Column tha contain the number of snp] [path to a folder
that contain your files] [options]

See string returned by @ref snp_charts.help_usage() for information

@subsection chart Charts
    - Quantitative bar chart (-q) : Show the number of gene (y) per number of snp (x)
    - Cumulative bar chart (-c): Show the number of gene (y) that have at least n snp (x)
    - Monoheatmap (-u): Cumulative bar chart but it's a heatmap
    - global heatmap (-g): concatenation of all Monoheatmap

@subsection format Formats
    - as tsv (-t)
    - as png (-k)
    - as svg (-b)
    - as a pop-up (-d)


@section authors Author(s)
  - Created by Marchal Florent on 06/05/2024.
This module has been made in 2024 during an internship at the UMR Agap, GE2pop (France, Montpellier)

@section libs Librairies/Modules
  - os
  - json
  - matplotlib.pyplot
  - sys
  - getopts_parser
  - utilities
"""
# Native libraries
import os
import json

try:
    from scripts.utilities.utilities import export_list_in_tsv_as_rows, chart_export, make_bar_char, make_heatmap, \
        parse_line, extract_data_from_table, FilterError
    from scripts.getopts_parser import getopts_parser


except ModuleNotFoundError as E:
    print(f"Module not found : {E}\n")
    if E.name == "getopts_parser":
        print(f"Make sure that 'getopts_parser/' is in the same folder as {__file__}")

    elif E.name == "utilities":
        print(f"Make sure that 'utilities/' is in the same folder as {__file__}")

    exit(1)

__author__ = "Marchal Florent"
__credits__ = ["Florent Marchal", "Vetea Jacot", "Concetta Burgarella", "Vincent Ranwez", "Nathalie Chantret"]


# List of arguments used to start this program with Their shorten name, (default value and type)
__getopts__ = {
    # long name:                (shorten_name, (default_value, type))
    #                                           None = (None, str)

    # Arguments:
    "name_column=":             ("n:", None),
    "snp_column=":              ("s:", None),
    "path=":                    ("p:", ("data", str)),
    "help": ("h", None),

    # Parameters :
    "file_separator=":          ("f:", ("\t", str)),
    "output_path=":             ("o:", ("output/", str)),
    "job_name=":                ("j:", ("Unnamed", str)),
    "max_length=":              ("m:", (20, int)),
    "show_values=":             ("e:", (None, int)),

    # long name:                (shorten_name, (inactive value, active value))
    #                                           None = (True, False)

    "output_warning":           ("w", (True, False)),
    "sort_by_name":             ("r", (True, False)),
    "simplified":               ("i", None),
    "global_heatmap":           ("g", None),
    "quantitative_barchart":    ("q", None),
    "cumulative_barchart":      ("c", None),
    "cumulative_heatmap":       ("u", None),
    "tsv":                      ("t", None),
    "png":                      ("k", None),
    "show":                     ("d", None),
    "svg":                      ("v", None),
    ("uniform_y", "uniform"):   ("y", None)

}


def help_usage():
    option_list = ""
    for key, (short_key, _) in __getopts__.items():
        if short_key[-1] != ":":
            option_list += f'{short_key} / {key}, \t'
        else:
            option_list += f'{short_key[:-1]} / {key}, \t'

    return ("python3 snp_analyser.py [Gene name Column] [Snp column] [Path to your files] [Options]"
            f"\nOptions are : "
            f"\n{option_list}"
            f"\nSee README.md for more details\n")


def filter_integer_greater_or_equal_to_0_ignore_0(key: str or int, dictionary: dict = None) -> bool:
    """!
    @brief Test if the value in front of the key @p key inside @p dictionary can be cast into an integer bigger than 0.

    Meant to be used inside  @ref extract_data_from_table as a "filter_" using a lambda.
        - value == 0 : Line is ignored
        - value > 0 : Line is kept
        - Value < 0 : Error is raised

    @warning does not support "1.0 nor "1,0" format

    @param key : str => A key contained by @p dictionary.
    @param dictionary : dict = None => A dictionary that contain @p key. If None, the function will assume that
    @code dictionary = {key: key} @endcode

    @return bool => True : Yes
                    False : value equal to 0
    @raise ValueError when the value is an integer lower than 0

    """
    if dictionary is None:
        dictionary = {key: key}

    value = int(dictionary[key])
    if value == 0:
        return False
    elif value < 0:
        raise ValueError(f"Integer lower or equal to 0 : {value}")
    return True


def compile_gene_snp(genes_snp: dict[str, any], dict_of_number: dict[int, dict[str, int]] = None,
                     group: str = "None") -> dict[int, dict[str, int]]:
    """!
    @brief Extract the number of snp of all genes contained in @p genes_snp (snp = @p genes_snp 's values).
    Each number of snp is stored inside a new dictionary (@p dict_of_number 's keys). A dict is created in front
    of all keys (i.e. snp number). This dict contain the @p group (key) and the number of occurrences of this
    snp number for this key.

    @param genes_snp : dict[str,any] => A dictionary from @ref extract_data_from_table.
        e.g. {gene_1: number_of_snp_in_gene_1} =>  @code {"gene1": 3}  @endcode
        @note Values (number of snp) inside this dict are trans typed into integers.
    @param dict_of_number : dict[int, dict[str,int]] = None.
        A dict with the same structure as dictionaries returned by this function.
    @param group : str = "None" => Each occurrence of a number of snp increment the counter related to this group.

    @return dict[int, dict[str, int]] => A dictionary that store all number of snp found along with the number of
    occurrences {number_of_snp_1 : {group1: number_of_occurrences_of_number_of_snp_1_in_this_group}

    @warning values @p genes_snp are cast into integer. Also, there is no verification made to see if the values are
    positive. We assume that data has been filtered using @ref filter_integer_greater_or_equal_to_0_ignore_0 in
    @ref extract_data_from_table
    """
    dict_of_number = {} if dict_of_number is None else dict_of_number

    for _, snp_count in genes_snp.items():
        snp_count = int(snp_count)

        # Add this 'snp_count' to dict_of_number
        if snp_count not in dict_of_number:
            dict_of_number[snp_count] = {}

        # Add this 'group' to dict_of_number[snp_count]
        if group not in dict_of_number[snp_count]:
            dict_of_number[snp_count][group] = 0

        dict_of_number[snp_count][group] += 1

    return dict_of_number


def make_data_matrix(compiled_dict : dict[int, dict[str, int]], group: str, *groups: str,
                     simplified: bool = True, max_length: int = None) -> (list[list[int]], list[int]):
    """!
    @brief Use a compiled dict from @ref compile_gene_snp to create a matrix of value.
    Each lines represent one group (@p groups & @p group).
    Each line contain the number of genes that contain at least n genes

    @parblock
    @note Only groups specified in @p groups and @p group are extracted from @p compiled_dict
    @endparblock

    @parblock
    @note For below example we will consider that :
        - @code compiled_dict = {1: {"E. coli": 5},
                           2: {"E. coli": 3, "HIV": 4},
                           3: {"HIV": 2}
                           4: {"E. coli": 1, "HIV": 1}
                           } @endcode
    @endparblock

    @param compiled_dict : dict[int,dict[str,int]] => a compiled dict from @ref compile_gene_snp
    @param group : str => A group name (e.g. Species name : "E. coli")
    @param *groups : str => Same as @p group. Additional species name.
    @param simplified : bool = True => Do number of snp represented by 0 gene are deleted from the result ?
    @note This the presence / absence of a number is affected by other groups
        - with group="E. coli":
            - If True: @code ([[5, 3, 1]], [1, 2, 4]) @endcode
            - If False: @code ([[5, 3, 0, 1]], [1, 2, 3, 4]) @endcode
        - with group="E. coli" and groups= ("HIV", ):
            - If True: @code ([[5, 3, 0, 1], [0, 4, 2, 1]], [1, 2, 4])  @endcode
    @param max_length : int = None => Limit the length of each lines of the matrix. should be greater than 0. If not,
    max_length is ignored.

    @return (list[list[int]], list[int]) => Return a matrix and a list. Each lines of the matrix represent the
    number of genes of a species (group) that contain n snp : If the position 1 of a line is equal to 4, there is
    4 genes in the selected species that contain list[i] snp.
    The y axes can be labeled using @p group and @p groups (order in conserved) and  the x-axis is labeled by
    the list[int].

    @note For some return example, see @p simplified for return example
    """
    if max_length is None or max_length <= 0:
        max_length = None

    # Variables
    groups = [group, *groups]                       # Merge @p group and @p groups
    data = [list() for _ in range(0, len(groups))]  # Data matrix
    last_x_value = 0                                # Remember the last value
    sorted_x_values = sorted(compiled_dict)         # Snp size from compiled_dict sorted y size

    # Apply the @p max_length by reducing the size of "sorted_x_values"
    if max_length is not None and len(sorted_x_values) >= max_length:
        sorted_x_values = sorted_x_values[:max_length]

    # Fill data
    for x_values in sorted_x_values:
        group_at_this_position = compiled_dict[x_values]

        # missing_x_values is append at by each groups in order to extend the matrix.
        # the last value (pos -1) is always rewrote in the next loop
        missing_x_values = [0]
        if simplified is False:
            # When there is snp numbers represented by 0 genes, we add them using missing_x_values
            missing_x_values *= (x_values - last_x_value)

        # For each group, extend the matrix
        for i, group_name in enumerate(groups):
            data[i].extend(missing_x_values)

            # Store the value of this snp number
            if group_name in group_at_this_position:
                data[i][-1] = group_at_this_position[group_name]

        # Update last_x_value
        last_x_value = x_values

    # Return the matrix and the legend
    if simplified is True:
        return data, sorted_x_values

    else:
        return data, list(range(1, sorted_x_values[-1] + 1))


def generate_cumulative_list(list_of_numbers: list[int] or list[float], reversed_=False) -> list[int]:
    """!
    @brief Take a list of number and sum all values.
    @code
    [0, 5, 6, 1] => [0+5+6+1, 5+6+1, 6+1, 1] == [12, 12, 7, 1]
    @endcode

    @param list_of_numbers : list[int] or list[float] => A list that contain numbers.
    @param reversed_ = False => Do the accumulation start at the end and end at the beginning.

    @return list[int] => A list of number

    """

    cumulative_list = []
    tot = 0

    # Select the correct way to browse list_of_numbers
    if reversed_:
        index_range = list(range(-1, -(len(list_of_numbers) + 1), -1))
    else:
        index_range = range(0, len(list_of_numbers))

    # Fill cumulative_list
    for i in index_range:
        tot += list_of_numbers[i]
        cumulative_list.append(tot)

    # Reverse cumulative_list when needed
    if reversed_:
        cumulative_list.reverse()

    return cumulative_list


def main(path: str, name_column: str, snp_column: str, file_separator: str = "\t",
         simplified: bool = True, max_length: int = None,
         output_path: str = "output", output_warning: bool = True, job_name: str = None,
         global_heatmap: bool = True, quantitative_barchart: bool = False, cumulative_barchart: bool = False,
         cumulative_heatmap: bool = False,
         tsv: bool = False, png: bool = False, show: bool = False, svg: bool = True,
         sort_by_name: bool = True, uniform_y: bool = True,
         show_values: int = -1) -> int:
    """!
    @brief Create a number of chart related to snp analysis.

    As an example we will consider the following flatFile :
    | GeneName | GeneID | NumberOfSnp | GeneSize |
    |----------|--------|-------------|----------|
    | Gene1    | 123    | 5           | 1000     |
    | Gene2    | 456    | 10          | 2000     |

    @param path : str => Path that lead to a number of flatfile :
        - .json : {complete file path : Species name}
        - folder : Use file inside the folder (does not scan the folder recursively)

    @param name_column : str => Name of the column that contain a primary key e.g. GeneName, GeneID. If two line
    have the same "primary key", the last one will be used.
    @param snp_column : str => Name of the column that contain a count of snp e.g. NumberOfSnp
    @param file_separator : str = "\t" => The separator used in all flat file considered.
    @param simplified : bool = True =>  Do number of snp represented by 0 gene are deleted from the result
    @param max_length : int = None => Maximum length of each graph. keep the nth first result. Should be greater or
    equal to 1 otherwise, it would be ignored.
    @param output_path : str = "output" => Where graphs are saved.
    @param output_warning : bool = True => Ask confirmation when at least one file can be erased by this program.
    @param job_name : str = None => A name for this execution. (This creates a separated folder in
     @p output_path). Default = "unnamed"
    @param global_heatmap : bool = True => Do a heatmap that is the combination of all cumulative_heatmap is created
    @param quantitative_barchart : bool = False => Do this program create a barchart of snp distribution for each
    file (Number of gene that have n snp)
    @param cumulative_barchart : bool = False => Do this program create a barchart of snp distribution for each
    file ? (Number of gene that have AT LEAST n snp)
    @param cumulative_heatmap : bool = False => Do this program create a heatmap of snp distribution for each
    file ? (Number of gene that have AT LEAST n snp)
    @param tsv : bool = False => Do values used for chart are saved in a flatfile (.tsv)
    @param png : bool = False => Do created charts are saved as png
    @param show : bool = False => Do created charts are saved are shown
        @warning Each time a chart is shown, the program stop. It will resume when the chart is closed.
    @param svg : bool = True => Do created charts are saved as svg (vectorize image)
    @param sort_by_name : bool = True => Do species are sorted in lexicographic order ?
    @param uniform_y : bool = True => Do all barchart share the same y axis ?
    @param show_values : int = None => If greater or equal to 0, all cells will contain theirs values. if lower than 0,
    text in cell in automatically determined (can be ugly when show is True, but assure that the text is good in
    png and svg). If None, nothing happen.

    @return int => if greater than 0, an error occurred.
    - 1 job stopped by user
    - 2 no species found
    """
    parameters = locals().copy()

    # Set a default name for output_path
    if job_name is None or len(job_name) == 0:
        job_name = "job_name"

    # Assure that max_length is None or greater or equal to 1
    if max_length <= 0:
        max_length = None
    if show_values <= 0:
        show_values = -1

    # Assure that at least one thing will be generated
    if not (global_heatmap or quantitative_barchart or cumulative_barchart or cumulative_heatmap):
        global_heatmap = True
    # Assure that at least one thing will be generated
    if not (tsv or svg or png or show):
        png = True

    # ---- ---- Path Management ---- ---- #
    # Assure that @p output_path point to a folder
    if output_path is None or output_path == "":
        output_path = "output/"

    if output_path[-1] not in ("/", "\\"):
        output_path += "/"

    # Create @p output_path if needed
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Load files
    if path[-5:] == ".json":
        with open(path, "r", encoding="utf-8") as js_flux:
            path_translation = dict(json.load(js_flux))
            list_of_files = path_translation
            file_path_prefix = ""

    else:
        # Assure that @p path point to a folder
        if path[-1] not in ("/", "\\"):
            file_path_prefix = path + "/"
        else:
            file_path_prefix = path

        list_of_files = os.listdir(file_path_prefix)
        path_translation = {}

    # Create file and directory path
    output_dir = output_path + job_name + "/"
    file_prefix = output_dir + job_name + "_"
    heatmap_prefix = file_prefix + "Heatmap"
    cumulative_prefix = file_prefix + "CumulativeBarchart_"
    quantitative_prefix = file_prefix + "QuantitativeBarchart_"

    # Generate output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Verify that the folder is empty
    elif output_warning and os.listdir(output_dir):
        r_ = input(f"Folder is not empty. Some files can be lost. ({output_dir})\nContinue ? (y / n) :")

        if r_.lower() not in ("y", "ye", "yes", "t", "tr", "tru", "true"):
            print("Job stopped")
            return 2

    # ---- ---- Export Control ---- ---  "
    heat_png = f"{heatmap_prefix}" if png and global_heatmap else None
    heat_tsv = f"{heatmap_prefix}" if tsv and global_heatmap else None
    heat_svg = f"{heatmap_prefix}" if svg and global_heatmap else None
    heat_show = show and global_heatmap

    c_heat_png = f"{heatmap_prefix}" if png and cumulative_heatmap else None
    c_heat_tsv = f"{heatmap_prefix}" if tsv and cumulative_heatmap else None
    c_heat_svg = f"{heatmap_prefix}" if svg and cumulative_heatmap else None
    c_heat_show = show and cumulative_heatmap

    c_bar_png = f"{cumulative_prefix}" if png and cumulative_barchart else None
    c_bar_tsv = f"{cumulative_prefix}" if tsv and cumulative_barchart else None
    c_bar_svg = f"{cumulative_prefix}" if svg and cumulative_barchart else None
    c_bar_show = show and cumulative_barchart

    q_bar_png = f"{quantitative_prefix}" if png and quantitative_barchart else None
    q_bar_tsv = f"{quantitative_prefix}" if tsv and quantitative_barchart else None
    q_bar_svg = f"{quantitative_prefix}" if svg and quantitative_barchart else None
    q_bar_show = show and quantitative_barchart

    # ---- ---- Load files ---- ----
    all_snp = {}                        # {Number_of_snp, {File_name : Number_of_genes_with_this_number_of_snp}
    all_species = []                    # List all targeted files

    # process all files and load snp into all_species
    for files in list_of_files:
        if files[0] == ".":
            continue

        try:
            files_dict = extract_data_from_table(f"{file_path_prefix}{files}", key=name_column, value=snp_column,
                                                 filter_=lambda _, val, dict_:
                                                            filter_integer_greater_or_equal_to_0_ignore_0(val, dict_),
                                                 # "_" in lambda is mendatory due to how extract_data_from_table works.
                                                 separator=file_separator)
        except FilterError as E:
            print(f"An error occurred : {E}\n"
                  
                  f"This program only accept positive integer in the following format : '1000', '1_000', '+1000'.\n"
                  f"Look for miss formated data in {snp_column}")

            return 5

        if files in path_translation:
            all_species.append(path_translation[files])
        else:
            all_species.append(files)

        all_snp = compile_gene_snp(files_dict, all_snp, group=all_species[-1])

    if sort_by_name:
        all_species.sort()

    # ---- ---- Matrix and chart generation ---- ----
    if not all_species:
        return 4

    data, x_legend = make_data_matrix(all_snp, *all_species, simplified=simplified, max_length=max_length)

    # Uniformize all y axis
    if uniform_y:
        max_quantitative_value = 0
        max_cumulative_value = 0
        for lines in data:
            max_quantitative_value = max(max(lines), max_quantitative_value)
            max_cumulative_value = max(sum(lines), max_quantitative_value)
        max_quantitative_value += 1
        max_cumulative_value += 1
    else:
        max_quantitative_value = None
        max_cumulative_value = None

    if len(x_legend) == x_legend[-1]:
        # if the legend is equivalent of the automatic one, we use the automatic legend
        # (e.g. when @p simplified is False or when there is no simplification),
        x_legend = None

    else:
        # When x_legend is not composed of str and @p simplified is True, BarChart have weird behaviour.
        x_legend = [str(item) for item in x_legend]

    for i in range(0, len(data)):
        line_name = all_species[i]

        # Make quantitative barchart
        if quantitative_barchart:
            make_bar_char(data[i], x_legend=x_legend, chart_name=line_name,
                          ylabel="Number of genes",
                          xlabel="Number of snp",
                          title=f"Number of snp per genes in {line_name}",
                          show=q_bar_show,
                          png=f"{q_bar_png}{line_name}" if q_bar_png else None,
                          tsv=f"{q_bar_tsv}{line_name}" if q_bar_tsv else None,
                          svg=f"{q_bar_svg}{line_name}" if q_bar_svg else None,
                          y_max_value = max_quantitative_value
                          )

        # Replace the quantitative list by a cumulative list
        data[i] = generate_cumulative_list(data[i], reversed_=True)

        # Make cumulative Barchart
        if cumulative_barchart:
            make_bar_char(data[i],
                          show=c_bar_show, x_legend=x_legend, chart_name=line_name,
                          ylabel="Number of genes",
                          xlabel="Number of snp",
                          title=f"Number of genes with at least n snp in {line_name}",
                          png=f"{c_bar_png}{line_name}" if c_bar_png else None,
                          tsv=f"{c_bar_tsv}{line_name}" if c_bar_tsv else None,
                          svg=f"{c_bar_svg}{line_name}" if c_bar_svg else None,
                          y_max_value = max_cumulative_value)

    # Heatmap generation
    if cumulative_heatmap:
        for i, lines in enumerate(data):
            make_heatmap([lines], y_legend=[all_species[i]], x_legend=x_legend,
                         title=f"Number of genes with at least n SNP : {all_species[i]}",
                         xlabel="Number of snp",
                         ylabel="Species names",
                         show=c_heat_show,
                         png=f"{c_heat_png}_{all_species[i]}" if c_heat_png else None,
                         tsv=f"{c_heat_tsv}_{all_species[i]}" if c_heat_tsv else None,
                         svg=f"{c_heat_svg}_{all_species[i]}" if c_heat_svg else None,
                         contain_number=show_values)

    if global_heatmap:
        make_heatmap(data, y_legend=all_species, x_legend=x_legend,
                     title="Number of genes with at least n SNP",
                     xlabel="Number of snp",
                     ylabel="Species names",
                     show=heat_show,
                     png=heat_png + "_global" if heat_png else None,
                     tsv=heat_tsv + "_global" if heat_tsv else None,
                     svg=heat_svg + "_global" if heat_svg else None,
                     contain_number=show_values)

    return 0


def main_using_getopts(argv: list[str] or str):
    """@brief start @ref snp_charts.main using a string or sys.argv[1:]
    Example :
        - main_using_getopts(sys.argv[1:])
        - main_using_getopts("name_column snp_column tests/TargetedFiles.json -m 20 -gv -w -j Tests -e -1")

    @param argv : list[str] or str => List of argument (strings). Usually sys.argv[1:].
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
    @return int => exit code
    """

    if isinstance(argv, str):
        argv = argv.split(" ")

    try:
        main_params = getopts_parser.getopts(argv, __getopts__, "name_column", "snp_column",
                                             help_options=("help", "h"),
                                             raise_errors=False,
                                             help_message=help_usage())
        if isinstance(main_params, int):
            exit(main_params + 10)

    except Exception as E:
        print(E)
        print(help_usage())
        return 5

    exit_code = main(**main_params)
    if exit_code == 2:
        print("Not enough species.")

    return exit_code


if __name__ == "__main__":
    import sys
    exit(main_using_getopts(sys.argv[1:]))



