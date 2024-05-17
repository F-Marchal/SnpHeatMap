try:
    import matplotlib.pyplot as plt

except ModuleNotFoundError as E:
    print(f"Module not found : {E}\n"
          f"Open a terminal and try : "
          f"\n\tpip install matplotlib"
          
          f"\n\nIf pip is not found, you can install it using : "
          f"\nOn linux or MacOs: "
          f"\n\tpython -m ensurepip --upgrade"
          f"\nOn Windows : "
          f"\n\tpy -m ensurepip --upgrade"
          )
    exit(1)


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


def chart_export(data: list[list[int]], show: bool = False, png: str = None, tsv: str = None, svg: str = None,
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
                  chart_name: str = None,
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
    @param chart_name : str = None => A name that will be used if tsv is not None to name a line.
    @param title : str = None => A title for this chart
    @param xlabel : str = None => A title for the x-axis
    @param ylabel : str = None => A title for the y-axis
    @param show : bool = False => Do current plot will be displayed ?
    @param png : str = None => Give a path to export the current plot as png
    @param tsv : str = None => Give a path to export @p data into a tsv.
    @param svg : str = None => Give a path to export the current plot as svg
    @param erase_last_plt : bool = True => If True, last plot is removed from @ref matplotlib.pyplot display
    """

    # Clear last plot
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

    # Assure that y-axis and x-axis use display only integer. (Just some plt magic)
    if y_legend_is_int:
        plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    if x_legend_is_int:
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Add labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.gcf().canvas.manager.set_window_title(title)

    # export chart
    chart_export(data=[data], x_legend=x_legend, tsv=tsv, png=png, show=show, svg=svg, y_legend=[chart_name])


def make_heatmap(data: list[list[int]],
                 x_legend: list = None, y_legend: list = None,
                 title: str = None, xlabel: str = None, ylabel: str = None,
                 show: bool = False, png: str = None, tsv: str = None, svg: str = None,
                 erase_last_plt: bool = True, contain_number: int = None,
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
    @param contain_number : int = None => If greater or equal to 0, all cells will contain theirs values. if lower than 0,
    text in cell in automatically determined. If None, nothing happen.
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
        plt.clf()
        plt.cla()
        plt.close('all')

    # Number of rows and columns
    num_rows = len(data)
    num_cols = len(data[0])

    # Create heatmap
    fig_size = (num_cols + 1, max(num_rows + 1, 4))

    plt.figure(figsize=fig_size)
    plt.imshow(data, cmap=cmap, interpolation='nearest', vmin=1)

    # Add ticks
    if x_legend:
        plt.xticks(range(1, num_cols+1), x_legend)
    else:
        x_legend = [str(i) for i in range(1, num_cols + 1)]
        plt.xticks(range(num_cols), x_legend)

    if y_legend:
        plt.yticks(range(num_rows), y_legend)

    # Add labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.gcf().canvas.manager.set_window_title(title)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Add color bar
    cbar = plt.colorbar()
    cbar.set_label('Number of genes')
    plt.gca().set_aspect('equal', adjustable='box')
    if contain_number is not None:
        # Add text labels inside heatmap cells
        units = ['', 'k', 'M', 'G', 'T', 'P']

        if contain_number < 0:
            font_size = min(12, min(fig_size) * 12)
        else:
            font_size = contain_number

        for i in range(len(data)):
            for j in range(len(data[0])):
                # pretreatment
                str_data = str(data[i][j])
                data_length = len(str_data)
                data_pos = data_length % 3
                data_units = data_length // 3

                # Round data at 3 significant numbers
                round_pos = max(data_length - 3, 0)
                round_pos *= -1
                str_data = str(round(data[i][j], round_pos))

                # Format units
                if data_units - 1 <= len(units):
                    # No unit
                    if data_units == 0:
                        str_data += "\n"    # Uniform text with other data

                    elif data_pos == 0:
                        str_data = f"{str_data[:3]}\n{units[data_units - 1]}"

                    else:
                        str_data = f"{str_data[:data_pos]},{str_data[data_pos: 3]}\n{units[data_units]}"

                else:
                    raise ValueError("Too many snp : " + str(data[i][j]))

                # Place text
                plt.text(j, i, f'{str_data}', ha='center', va='center', color=test_color, fontsize=font_size)

    chart_export(data=data, y_legend=y_legend, x_legend=x_legend, tsv=tsv, png=png, show=show, svg=svg)


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
                raise ValueError(f"Both key ('{key}') and value('{value}') should be contained inside the legend : {legend}")

            if legend[-1][-1] == "\n":
                legend[-1] = legend[-1][:-1]
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
