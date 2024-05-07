
import os
import json
import matplotlib.pyplot as plt
import math

def parse_line(legend: list[str], line: str, separator: str = "\t") -> dict[str, str]:
    """! @brief Turn a line form a flat File with its legend and turn it into a dictionary.
    @param legend : Names all the line's columns. Example: ["A", "B", "C"]
    @param line : Contains all the line's values. Example: "1|2|3|", "1|2|3" ...
    @param separator : The symbol that splits the line's values. Example: "|", "\n", "\t" ...
    @return A dictionary composed of legend's values and line's values.
        Example:  @code {"A": "1", "B": "2", "C": "3"}  @endcode (using previous examples)

    @note The returned dict always contain the same umber of object than legend.
        - If legend > line : part of the legend's values will point to an empty string
        - If legend < line : part of the line will be ignored
    """

    parsed_line = {}
    split_line = line.split(separator)  # Split the line in "columns"

    # Fill parsed_line
    for i, column_name in enumerate(legend):
        cell_value = split_line[i] if i < len(split_line) else ""
        parsed_line[column_name] = cell_value
        i += 1

    return parsed_line


def extract_data_from_table(path: str, key: str, value: str, separator: str = "\t",
                            legend: list = None, filter_: callable = None) -> dict[str, any]:
    """!
    @brief Read a table contained inside a flatFile (e.g. tsv, csv, ...)

    @warning If the column @p key contains the same value multiple times, only the last one is kept.

    @param path : Path to a flatFile.
    @param key : A column name that can be found in the legend. This will be used as a key in the returned dict.
            Example: "Column3"
    @param value : A column name that can be found in the legend. This will be used as a value in the returned dict
            WHEN @p filter returns None or True. Example: "Column2"
    @param separator : The symbol that splits line's values. Example: "|", "\n", "\t" ...
    @param legend : If None: The first non-empty line in the file split using @p separator. Else: A list of column
            names. Example: [Column1, Column2, Column3]
    @param filter_ : A function that accepts 3 arguments: @p key, @p value, and the parsed line (dict). It
            selects/generates the value present next to each key.
                - If it returns True or None: value in the column @p value.
                - If it returns False: this line is ignored.
                - Else: The returned value is used (instead of the content of the column @p value).
    @note filter_ is called one time per line.
    @return A dictionary: {values in the column @p key (values that do not pass @p filter_ are ignored): values in the column @p value OR value returned by @p filter_}
    """
    # Open file
    flux = open(path, "r", encoding="UTF-8")
    data = {}

    # Researched keys and values should be contained inside the legend.
    if legend is not None and (key not in legend or value not in legend):
        raise ValueError(f"Both key ('{key}') and value {value} should be contained inside legend : {legend}")

    # Fill data
    for line in flux:
        # Special cases : empty line | empty file
        if line == "\n" or line == "":
            continue

        # Special cases : Unknown legend
        if legend is None:
            legend = line.split(separator)

            # Researched keys and values should be contained inside the legend.
            if key not in legend or value not in legend:
                raise ValueError(f"Both key ('{key}') and value {value} should be contained inside the legend : {legend}")

            continue
        
        # Parse the line
        parsed_line = parse_line(legend, line, separator)
        
        # Apply the filter_
        func_result = filter_(key, value, parsed_line) if filter_ else None
        if func_result is False:
            continue

        if func_result is None or func_result is True:
            # Save  parsed_line[value]
            data[parsed_line[key]] = parsed_line[value]

        else:
            # Save  func_result
            data[parsed_line[key]] = func_result

    # Close file
    flux.close()

    return data


def greater_than_0_int_filter(_, key: int=None, dictionary: dict = None) -> bool:
    """!
    @brief Test if the value in front of the key @p key inside @p dictionary can be an integer bigger than 0.
    Meant to be used inside  @ref extract_data_from_table as a "filter_"

    @param _ => Unused parameter. Exist due to how "filter_" in @extract_data_from_table works
    @param key : int = None => A key contained by @p dictionary.
    @param dictionary : dict = None => A dictionary that contain @p key

    @return bool => True : Yes ; False : No.

    """
    try:
        return int(dictionary[key]) > 0
    except ValueError:
        return False


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
    @param max_length : int = None => Limit the length of each lines of the matrix.

    @return (list[list[int]], list[int]) => Return a matrix and a list. Each lines of the matrix represent the
    number of genes of a species (group) that contain n snp : If the position 1 of a line is equal to 4, there is
    4 genes in the selected species that contain list[i] snp.
    The y axes can be labeled using @p group and @p groups (order in conserved) and  the x-axis is labeled by
    the list[int].

    @note For some return example, see @p simplified for return example
    """
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
    @encode

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


def export_list_in_tsv_as_rows(path: str, *rows, file_mode="w", encoding="UTF-8",
                               y_legend: list = None, x_legend: list = None):
    """!
    @brief Accept a number of list that represent rows of a tab and turn it intoo a tsv (flat file).

    @warning  Any "\t" or "\n" in rows' values will disrupt lines and / or columns

    Parameters :
        @param path : str => path (and name) of the file that will be written.
        @param *rows => A number of list
        @param file_mode = "w" => "w" or "a"
            - "w": if a file with the same path exist, old file is erased
            - "a": if a file with the same path exist, old file append new values
        @param encoding = "UTF-8" => File encoding
        @param y_legend : list = None => A list of item to be display in the first column
        @param x_legend : list = None => A list of item to be display at top of the file

    """
    # WARNING: \t and \n in rows can distributed lines / columns

    # Open file
    file_flux = open(path, mode=file_mode, encoding=encoding)

    if x_legend:
        # Add the legend
        rows = [x_legend, *rows]

        # Assure that y_legend is aligned with the correct line
        # (x_legend create a new line that shift y_legend)
        if y_legend is not None and len(y_legend) < len(rows):
            y_legend = ["", *y_legend]

    i = 0   # Initialise i for the last block of instruction
    for i, lines in enumerate(rows):

        # Add y_legend at the beginning of each lines
        if y_legend:
            if i < len(y_legend):   # Assure that y_legend can not create errors
                file_flux.write(y_legend[i])

            file_flux.write("\t")

        # Write line content
        for word in lines:
            file_flux.write(str(word) + "\t")

        # End line
        file_flux.write("\n")

    # Assure that y_legend is completely written
    if y_legend:
        while i < len(y_legend) - 1:
            file_flux.write(y_legend[i])
            file_flux.write("\t")
            i += 1

    # Close file
    file_flux.close()


def _chart_export(data: list[list[int]], show: bool = False, png: str = None, tsv: str = None, svg: str = None,
                  x_legend: list = None, y_legend: list = None):
    """!
    @brief Export the current chart.

    @note if data is the only argument, nothing happen.

    @param data : list[list[int]] => A matrix of values
    @param show : bool = False => Do current plot will be displayed ?
    @param png : str = None => Give a path to export the current plot as png
    @param tsv : str = None => Give a path to export @p data into a tsv.
    @param svg : str = None => Give a path to export the current plot as svg
    @param y_legend : list = None => When tsv is not none: A list of item to be display in the first column
                (@ref export_list_in_tsv_as_rows)
    @param x_legend : list = None => When tsv is not none: A list of item to be display at top of the file
                (@ref export_list_in_tsv_as_rows)

    """
    # Png export
    if png is not None:
        plt.savefig(png + ".png", format='png')

    # svg export (Scalable Vector Graphic )
    if svg is not None:
        plt.savefig(svg + ".svg", format='svg')

    # Export tsv (flat file)
    if tsv is not None:
        export_list_in_tsv_as_rows(tsv + ".tsv", *data, y_legend=y_legend, x_legend=x_legend)

    # show chart
    if show:
        plt.show()


def make_bar_char(data: list[int],
                  x_legend: list = None, x_legend_is_int: bool = True, y_legend_is_int: bool = True,
                  title: str = None, xlabel: str = None, ylabel: str = None,
                  show: bool = False, png: str = None, tsv: str = None, svg: str = None,
                  erase_last_plt: bool = True):
    """!
    @brief Create a @ref plt.bar using a bunch of argument.
    This function is made to assure a correct looking legend when used for snp.

    @param data : list[int] => A list of integer
    @param x_legend : list = None => Values used to legend the x-axis.
    @param x_legend_is_int : bool = True => Do x-axis represent oly integer
    @param y_legend_is_int : bool = True => Do y-axis represent oly integer
    @param title : str = None => A title for this chart
    @param xlabel : str = None => A title for the x-axis
    @param ylabel : str = None => A title for the y-axis
    @param show : bool = False => Do current plot will be displayed ?
    @param png : str = None => Give a path to export the current plot as png
    @param tsv : str = None => Give a path to export @p data into a tsv.
    @param svg : str = None => Give a path to export the current plot as svg
    @param erase_last_plt : bool = True => If True, last plot is removed from @ref matplotlib.pyplot memory
    """

    # Clear the last plot
    if erase_last_plt:
        plt.close('all')
        plt.clf()
        plt.cla()

    # Add ticks
    if x_legend:
        plt.bar(x_legend, data, color='skyblue')
        plt.xticks(range(len(data)), x_legend)

    else:
        x_legend = list(range(1, len(data) + 1))
        plt.bar(x_legend, data, color='skyblue')

    if y_legend_is_int:
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    if x_legend_is_int:
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Add labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    # export chart
    _chart_export(data=[data], x_legend=x_legend, tsv=tsv, png=png, show=show, svg=svg)


def make_heatmap(data: list[list[int]],
                 x_legend: list = None, y_legend: list = None,
                 title: str = None, xlabel: str = None, ylabel: str = None,
                 show: bool = False, png: str = None, tsv: str = None, svg: str = None,
                 erase_last_plt: bool = True, contain_number: bool = True,
                 test_color: str = "#a0a0a0", cmap: str = "jet",
                 ):
    """!
    @brief Create a heatmap using a bunch of argument.
    This function is made to assure a correct looking legend when used for snp.

    @param data : list[int] => A list of integer
    @param x_legend : list = None => Values used to label the x-axis.
    @param y_legend : list = None => Values used to label the y-axis.
    @param title : str = None => A title for this chart
    @param xlabel : str = None => A title for the x-axis
    @param ylabel : str = None => A title for the y-axis
    @param show : bool = False => Do current plot will be displayed ?
    @param png : str = None => Give a path to export the current plot as png
    @param tsv : str = None => Give a path to export @p data into a tsv.
    @param svg : str = None => Give a path to export the current plot as svg
    @param erase_last_plt : bool = True => If True, last plot is removed from @ref matplotlib.pyplot memory
    @param contain_number : bool = True => If True, all cells will contain theirs values.
    @param test_color : str = #a0a0a0 => HTML color code for text inside cells
        @note Only when contain_number is True
    @param cmap : str = jet => Color mod. supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG',
    'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
    'Grays', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r',
    'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn',
    'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy',
    'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2',
     'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', '
     YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary',
     'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r',
     'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth',
     'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_grey', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',
     'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg',
     'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'hot', 'hot_r', 'hsv', 'hsv_r',
     'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean',
     'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic',
     'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b',
     'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r',
     'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'
    """

    # Clear the last plot
    if erase_last_plt:
        plt.close('all')
        plt.clf()
        plt.cla()

    # Number of rows and columns
    num_rows = len(data)
    num_cols = len(data[0])

    # Create heatmap
    plt.figure(figsize=(num_cols + 1, max(num_rows + 1, 4)))
    plt.imshow(data, cmap=cmap, interpolation='nearest', vmin=1)

    # Add ticks
    if x_legend:
        plt.xticks(range(1, num_cols+1), x_legend)
    else:
        x_legend = list(range(1, num_cols + 1))
        plt.xticks(range(num_cols), x_legend)

    if y_legend:
        plt.yticks(range(num_rows), y_legend)

    # Add labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Add color bar
    cbar = plt.colorbar()
    cbar.set_label('Number of genes')

    if contain_number:
        # Add text labels inside heatmap cells
        units = ['', 'k', 'M', 'G', 'T', 'P']
        for i in range(len(data)):
            for j in range(len(data[0])):
                # pretreatment
                str_data = str(data[i][j])
                data_length = len(str_data)
                data_pos = data_length % 3
                data_units = data_length // 3

                # Format units
                if data_units - 1 <= len(units):
                    if data_units == 0:
                        pass

                    elif data_pos == 0:
                        str_data = f"{str_data[:3]}\n{units[data_units - 1]}"

                    else:
                        str_data = f"{str_data[:data_pos]},{str_data[data_pos: 3]}\n{units[data_units]}"

                else:
                    raise ValueError("Too many snp : " + str(data[i][j]))

                # Place text
                plt.text(j, i, f'{str_data}', ha='center', va='center', color=test_color)

        pass

    _chart_export(data=data, y_legend=y_legend, x_legend=x_legend, tsv=tsv, png=png, show=show, svg=svg)


def main(path: str, name_column: str, snp_column: str, file_separator: str = "\t",
         simplified: bool = True, max_length: int = None,
         output_path: str = "output", output_warning: bool = True, job_name: str = "unnamed",
         global_heatmap: bool = True, quantitative_barchart: bool = False, cumulative_barchart: bool = False,
         cumulative_heatmap: bool = False,
         tsv: bool = False, png: bool = False, show: bool = False, svg: bool = True,
         sort_names: bool = True) -> int:
    # WARING: Non recursiv, pas de trie

    # ---- ---- Path Management ---- ---- #
    # Assure that @p output_path point to a folder
    if output_path[-1] not in ("/", "\\"):
        output_path += "/"

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

    # Create @p output_path if needed
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Create file and directory path
    output_dir = output_path + "/" + job_name + "/"
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
            return 1

    # ---- ---- Export Control ---- ---  "
    heat_png = f"{heatmap_prefix}" if png and global_heatmap else None
    heat_tsv = f"{heatmap_prefix}" if tsv and global_heatmap else None
    heat_svg = f"{heatmap_prefix}" if svg and global_heatmap else None
    heat_show = show and global_heatmap

    c_heat_png = f"{heatmap_prefix}" if png and global_heatmap else None
    c_heat_tsv = f"{heatmap_prefix}" if tsv and global_heatmap else None
    c_heat_svg = f"{heatmap_prefix}" if svg and global_heatmap else None
    c_heat_show = show and global_heatmap

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

    # process all files and load snp into all_files
    for files in list_of_files:
        files_dict = extract_data_from_table(f"{file_path_prefix}{files}", key=name_column, value=snp_column,
                                             filter_=greater_than_0_int_filter, separator=file_separator)

        if files in path_translation:
            all_species.append(path_translation[files])
        else:
            all_species.append(files)

        all_snp = compile_gene_snp(files_dict, all_snp, group=all_species[-1])

    if sort_names:
        all_species.sort()

    # ---- ---- Matrix and chart generation ---- ----
    data, x_legend = make_data_matrix(all_snp, *all_species, simplified=simplified, max_length=max_length)

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
            make_bar_char(data[i], x_legend=x_legend,
                          ylabel="Number of genes",
                          xlabel="Number of snp",
                          title=f"Number of snp per genes in {line_name}",
                          show=q_bar_show,
                          png=f"{q_bar_png}{line_name}" if q_bar_png else None,
                          tsv=f"{q_bar_tsv}{line_name}" if q_bar_tsv else None,
                          svg=f"{q_bar_svg}{line_name}" if q_bar_svg else None,
                          )

        # Replace the quantitative list by a cumulative list
        data[i] = generate_cumulative_list(data[i], reversed_=True)

        # Make cumulative Barchart
        if cumulative_barchart:
            make_bar_char(data[i],
                          show=c_bar_show, x_legend=x_legend,
                          ylabel="Number of genes",
                          xlabel="Number of snp",
                          title=f"Number of genes with at least n snp in {line_name}",
                          png=f"{c_bar_png}{line_name}" if c_bar_png else None,
                          tsv=f"{c_bar_tsv}{line_name}" if c_bar_tsv else None,
                          svg=f"{c_bar_svg}{line_name}" if c_bar_svg else None,
                          )

    # Heatmap generation
    if cumulative_heatmap:
        for i, lines in enumerate(data):
            make_heatmap([lines], y_legend=[""], x_legend=x_legend,
                         title=f"Number of genes with at least n SNP : {all_species[i]}",
                         xlabel="Number of snp",
                         ylabel="Species names",
                         show=c_heat_show,
                         png=f"{c_heat_png}_{all_species[i]}" if c_heat_png else None,
                         tsv=f"{c_heat_tsv}_{all_species[i]}" if c_heat_tsv else None,
                         svg=f"{c_heat_svg}_{all_species[i]}" if c_heat_svg else None,
                        )

    if global_heatmap:
        make_heatmap(data, y_legend=all_species, x_legend=x_legend,
                     title="Number of genes with at least n SNP",
                     xlabel="Number of snp",
                     ylabel="Species names",
                     show=heat_show,
                     png=heat_png + "_global" if heat_png else None,
                     tsv=heat_tsv + "_global" if heat_tsv else None,
                     svg=heat_svg + "_global" if heat_svg else None,
                     )

    return 0


if __name__ == "__main__":
    main("TargetedFiles.json",
         "Contig_name", "BiAllelic_SNP", output_warning=False,
         max_length=20, tsv=False, cumulative_heatmap=True,

         job_name="Data")

    #TODO: Unitary test