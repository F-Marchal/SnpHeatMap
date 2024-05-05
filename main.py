
import os
import sys
import matplotlib.pyplot as plt
import math


def parse_line(legend: list[str], line: str, separator: str = "\t") -> dict[str, str]:
    """!
    @brief Converts : a line from a table to a dictionary using a specified legend and separator.
    @param legend : Names all the line's columns. Example: ["A", "B", "C"]
    @param line : Contains all the line's values. Example: "1|2|3|", "1|2|3" ...
    @param separator : The symbol that splits the line's values. Example: "|", "\n", "\t" ...
    @return A dictionary composed of legend's values and line's values.
        Example: {"A": "1", "B": "2", "C": "3"} (using previous examples)
    """
    parsed_line = {}
    split_line = line.split(separator)

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

    Parameters:
        @param path : Path to a flatFile.
        @param key : A column name that can be found in the legend. This will be used as a key in the returned dict. Example: "Column3"
        @param value : A column name that can be found in the legend. This will be used as a value in the returned dict WHEN @p filter returns None or True. Example: "Column2"
        @param separator : The symbol that splits line's values. Example: "|", "\n", "\t" ...
        @param legend : If None: The first non-empty line in the file split using @p separator. Else: A list of column names. Example: [Column1, Column2, Column3]
        @param filter_ : A function that accepts 3 arguments: @p key, @p value, and the parsed line (dict). It selects/generates the value present next to each key.
            - If it returns True or None: value in the column @p value.
            - If it returns False: this line is ignored.
            - Else: The returned value is used (instead of the content of the column @p value).
            @note filter_ is called one time per line.
    @return A dictionary: {values in the column @p key (values that do not pass @p filter_ are ignored): values in the column @p value OR value returned by @p filter_}
    """

    flux = open(path, "r", encoding="UTF-8")
    data = {}
    
    if legend is not None and (key not in legend or value not in legend):
        raise ValueError(f"Both key ('{key}') and value {value} should be contained inside legend : {legend}")

    for line in flux:
        # Special cases : empty line | empty file
        if line == "\n" or line == "":
            continue

        # Special cases : Unknown legend
        if legend is None:
            legend = line.split("\t")

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
            data[parsed_line[key]] = parsed_line[value]

        else:
            data[parsed_line[key]] = func_result

    flux.close()
    return data


def greater_than_0_int_filter(_, key: int=None, dictionary: dict = None) -> bool:
    """!
    @brief Test if the value in front of the key @p key inside @p dictionary can be an integer bigger than 0.
    Meant to be used inside  @extract_data_from_table as a "filter_"

    Parameters : 
        @param _ => Unused parameter. Exist due to how "filter_" in @extract_data_from_table works
        @param key : int = None => A key contained by @p dictionary.
        @param dictionary : dict = None => A dictionary that contain @p key
    
    Return : 
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

    Parameters :
        @param genes_snp : dict[str,any] => A dictionary from @ref extract_data_from_table.
            e.g. {gene_1: number_of_snp_in_gene_1} => {"gene1": 3}
            @note Values (number of snp) inside this dict are trans typed into integers.
        @param dict_of_number : dict[int, dict[str,int]] = None.
            A dict with the same structure as dictionaries returned by this function.
        @param group : str = "None" => Each occurrence of a number of snp increment the counter related to this group.
    Returns :
        @return dict[int, dict[str, int]] => A dictionary that store all number of snp found along with the number of
        occurrences {number_of_snp_1 : {group1: number_of_occurrences_of_number_of_snp_1_in_this_group}
    """
    dict_of_number = {} if dict_of_number is None else dict_of_number

    for _, snp_count in genes_snp.items():
        snp_count = int(snp_count)

        if snp_count not in dict_of_number:
            dict_of_number[snp_count] = {}

        if group not in dict_of_number[snp_count]:
            dict_of_number[snp_count][group] = 0

        dict_of_number[snp_count][group] += 1

    return dict_of_number


def make_data_matrix(compiled_dict : dict[int, dict[str, int]], group : str, *groups : str,
                     simplified: bool = True, max_length: int=None) -> (list[list[int]], list[int]):
    groups = [group, *groups]
    data = [list() for _ in range(0, len(groups))]
    last_x_value = 0
    sorted_x_values = sorted(compiled_dict)

    if max_length is not None and len(sorted_x_values) >= max_length:
        sorted_x_values = sorted_x_values[:max_length]

    for x_values in sorted_x_values:
        group_at_this_position = compiled_dict[x_values]

        missing_x_values = [0]
        if simplified is False:
            missing_x_values *= (x_values - last_x_value)

        for i, group_name in enumerate(groups):
            data[i].extend(missing_x_values)

            print(group_name, data[i])
            if group_name in group_at_this_position:
                data[i][-1] = group_at_this_position[group_name]

        last_x_value = x_values

    if simplified is True:
        return data, sorted_x_values

    else:
        return data, list(range(1, sorted_x_values[-1] + 1))



def generate_cumulative_list(list_of_numbers: list[int] or list[float], reversed_=False) -> list[int]:
    cumulative_list = []
    tot = 0

    if reversed_:
        index_range = list(range(-1, -(len(list_of_numbers) + 1), -1))
    else:
        index_range = range(0, len(list_of_numbers))

    for i in index_range:
        tot += list_of_numbers[i]
        cumulative_list.append(tot)

    if reversed_:
        cumulative_list.reverse()

    return cumulative_list

def export_list_in_tsv_as_columns(path: str, *columns, file_mode="w", encoding="UTF-8",
                                  y_legend: list=None, x_legend: list=None):
    file_flux = open(path, mode=file_mode, encoding=encoding)
    if y_legend is None:
        columns = [y_legend, *columns]
    row_index = 0
    someone_was_not_empty = True

    if x_legend:
        for item in x_legend:
            file_flux.write(item + "\t")
    file_flux.write("\n")

    while someone_was_not_empty:
        someone_was_not_empty = False

        for current_columns in columns:
            if len(columns) > row_index:
                file_flux.write(current_columns[row_index])
                someone_was_not_empty = True
            file_flux.write("\t")

        file_flux.write("\n")
        row_index += 1

    file_flux.close()

def export_list_in_tsv_as_rows(path: str, *rows, file_mode="w", encoding="UTF-8",
                                  y_legend: list=None, x_legend: list=None):
    file_flux = open(path, mode=file_mode, encoding=encoding)

    if x_legend:
        rows = [x_legend, *rows]

    i = 0
    for i, lines in enumerate(rows):
        if y_legend:
            if i < len(y_legend):
                file_flux.write(y_legend[i])
            file_flux.write("\t")

        for word in lines:
            file_flux.write(str(word) + "\t")
        file_flux.write("\n")

    if y_legend:
        while i < len(y_legend):
            file_flux.write(y_legend[i])
            file_flux.write("\t")

    file_flux.close()

def _chart_export(data: list[list[int]], show: bool = False, png: str = None, tsv: str = None, svg: str = None):
    # chart display
    if show:
        plt.show()

    if png is not None:
        plt.savefig(png + ".png", format='png')

    if svg is not None:
        plt.savefig(svg + ".svg", format='svg')

    if tsv is not None:
        pass
        # export_list_in_tsv_as_rows(tsv + ".tsv", *data, values)


def _chart_labels(x_positions, legend, title: str = None, xlabel: str = None, ylabel: str = None):

    plt.xticks(x_positions, legend)
    plt.xlabel(xlabel)
    plt.xlabel(ylabel)
    plt.title(title)


def make_bar_chart(values: list[int or float], legend: list[any],
                   title: str = None, xlabel: str = None, ylabel: str = None,
                   show: bool = False, png: str = None, tsv: str = None, erase: bool = True) -> plt:

    if erase is True:
        plt.clf()
        plt.close()

    # Create the bar chart
    x_positions = range(len(legend))
    plt.plot(x_positions, values, color='blue', alpha=0.7)

    # Set labels
    _chart_labels(x_positions, legend, xlabel=xlabel, ylabel=ylabel, title=title)
    _chart_export(legend=legend, values=values, tsv=tsv, png=png, show=show)

    return plt

def make_heatmap(data: list[list[int]], x_legend: list=None, y_legend: list=None):
    # Number of rows and columns
    num_rows = len(data)
    num_cols = len(data[0])

    # Create heatmap
    plt.imshow(data, cmap='viridis', interpolation='nearest')

    # Add ticks and labels
    if x_legend:
        plt.xticks(range(num_cols), x_legend)

    if y_legend:
        plt.yticks(range(num_rows), y_legend)

    # Add color bar
    plt.colorbar()

    # Show plot
    plt.show()
def main(path: str):
    all_snp = {}
    all_files = os.listdir(path)
    for files in all_files:
        files_dict = extract_data_from_table(path + "/" + files, key="Contig_name", value="BiAllelic_SNP",
                                             filter_=greater_than_0_int_filter)
        all_snp = compile_gene_snp(files_dict, all_snp, group=files)

    data, x_legend = make_data_matrix(all_snp, *all_files, simplified=True, max_length=28)
    for i in range(0, len(data)):
        data[i] = generate_cumulative_list(data[i], reversed_=True)

    if len(x_legend) == x_legend[-1]:
        x_legend = None
    make_heatmap(data, y_legend=all_files, x_legend=x_legend)


if __name__ == "__main__":
    main("tests")
    
    #TODO: Save function
    #TODO: Barplot for each files
    #TODO: Documentation
    #TODO: Test unitaire