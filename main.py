
import os
import json
import matplotlib.pyplot as plt


def parse_line(legend: list[str], line: str, separator: str = "\t") -> dict[str, str]:
    """! @brief Turn a line form a flat File with its legend and turn it into a dictionary.

    @cond Parameters:                                                                                       @endcond
        @param legend : Names all the line's columns. Example: ["A", "B", "C"]
        @param line : Contains all the line's values. Example: "1|2|3|", "1|2|3" ...
        @param separator : The symbol that splits the line's values. Example: "|", "\n", "\t" ...

    @cond Returns:                                                                                          @endcond
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

    @cond Parameters:                                                                                       @endcond
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

    @cond Returns:                                                                                          @endcond
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

    @cond Parameters:                                                                                       @endcond
        @param _ => Unused parameter. Exist due to how "filter_" in @extract_data_from_table works
        @param key : int = None => A key contained by @p dictionary.
        @param dictionary : dict = None => A dictionary that contain @p key

    @cond Returns:                                                                                          @endcond
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

    @cond Parameters:                                                                                       @endcond
        @param genes_snp : dict[str,any] => A dictionary from @ref extract_data_from_table.
            e.g. {gene_1: number_of_snp_in_gene_1} =>  @code {"gene1": 3}  @endcode
            @note Values (number of snp) inside this dict are trans typed into integers.
        @param dict_of_number : dict[int, dict[str,int]] = None.
            A dict with the same structure as dictionaries returned by this function.
        @param group : str = "None" => Each occurrence of a number of snp increment the counter related to this group.

    @cond Returns:                                                                                          @endcond
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

    @cond Parameters:                                                                                       @endcond
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

    @cond Returns:                                                                                          @endcond
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
    @brief [Function's description]


    Parameters :
        @param list_of_numbers : list[int]orlist[float] => [description]
        @param reversed_ = False => [description]

    Returns :
        @return list[int] => [description]

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
        plt.bar(range(1, len(data) + 1), data, color='skyblue')

    if y_legend_is_int:
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    if x_legend_is_int:
        # plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        pass

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
                 ):
    # Clear the last plot
    if erase_last_plt:
        plt.close('all')
        plt.clf()
        plt.cla()

    # Number of rows and columns
    num_rows = len(data)
    num_cols = len(data[0])

    # Create heatmap
    plt.imshow(data, cmap='viridis', interpolation='nearest', vmin=1)

    # Add ticks
    if x_legend:
        plt.xticks(range(1, num_cols+1), x_legend)
    else:
        plt.xticks(range(num_cols), range(1, num_cols + 1))

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
        for i in range(len(data)):
            for j in range(len(data[0])):
                plt.text(j, i, str(data[i][j]), ha='center', va='center', color='black')

    _chart_export(data=data, y_legend=y_legend, x_legend=x_legend, tsv=tsv, png=png, show=show, svg=svg)


def main(path: str, name_column: str, snp_column: str,
         simplified: bool = True, max_length: int = None,
         output_path: str = "output", output_warning: bool = True, job_name: str = "unnamed",
         heatmap: bool = True, barchart: bool = True,
         tsv: bool = False, png: bool = False, show: bool = False, svg: bool = True,
         file_name_to_species_path: str = None) -> int:
    # WARING: Non recursiv, pas de trie

    # ---- ---- Path Management ---- ---- #
    # Assure that @p output_path point to a folder
    if output_path[-1] not in ("/", "\\"):
        output_path += "/"

    # Assure that @p path point to a folder
    if path[-1] not in ("/", "\\"):
        path += "/"

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

    # File name to species name translation dict
    translation_dict = {}
    if file_name_to_species_path:
        if os.path.isfile(file_name_to_species_path):
            with open(file_name_to_species_path, "r") as file_flux:
                translation_dict = dict(json.load(file_flux))

        else:
            print(f"Unable to load {file_name_to_species_path}. File names will not be translated to species name.")

    # Verify that the folder is empty
    elif output_warning and os.listdir(output_dir):
        r_ = input(f"Folder is not empty. Some files can be lost. ({output_dir})\nContinue ? (y / n) :")

        if r_.lower() not in ("y", "ye", "yes", "t", "tr", "tru", "true"):
            print("Job stopped")
            return 1

    # ---- ---- Export Control ---- ---  "
    heat_png = f"{heatmap_prefix}" if png and heatmap else None
    heat_tsv = f"{heatmap_prefix}" if tsv and heatmap else None
    heat_svg = f"{heatmap_prefix}" if svg and heatmap else None
    heat_show = show and heatmap

    c_bar_png = f"{cumulative_prefix}" if png and barchart else None
    c_bar_tsv = f"{cumulative_prefix}" if tsv and barchart else None
    c_bar_svg = f"{cumulative_prefix}" if svg and barchart else None
    c_bar_show = show and barchart

    q_bar_png = f"{quantitative_prefix}" if png and barchart else None
    q_bar_tsv = f"{quantitative_prefix}" if tsv and barchart else None
    q_bar_svg = f"{quantitative_prefix}" if svg and barchart else None
    q_bar_show = show and barchart

    # ---- ---- Load files ---- ----
    all_snp = {}                        # {Number_of_snp, {File_name : Number_of_genes_with_this_number_of_snp}
    all_species = []                    # List all targeted files

    # process all files and load snp into all_files
    for files in os.listdir(path):
        files_dict = extract_data_from_table(f"{path}{files}", key=name_column, value=snp_column,
                                             filter_=greater_than_0_int_filter)

        if files in translation_dict:
            all_species.append(translation_dict[files])
        else:
            all_species.append(files)

        all_snp = compile_gene_snp(files_dict, all_snp, group=all_species[-1])

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
        if barchart:
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
        if barchart:
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
    if heatmap:
        make_heatmap(data, y_legend=all_species, x_legend=x_legend,
                     title="Number of genes with at least n SNP",
                     xlabel="Number of snp",
                     ylabel="Species names",
                     show=heat_show,
                     png=heat_png,
                     tsv=heat_tsv,
                     svg=heat_svg,
                     )

    return 0
if __name__ == "__main__":
    main("tests/data", "Contig_name", "BiAllelic_SNP", output_warning=False,
         max_length=5, file_name_to_species_path="tests/tranlateFileName.json", job_name="Example")




    #TODO: Documentation
    #TODO: Test unitaire