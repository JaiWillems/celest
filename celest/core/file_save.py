

import numpy as np


NEWLINE_STRING = "\n"
BLANK_LINE_STRING = "\n\n"
COLUMN_DELIMITER = "\t"


def _save_data_as_txt(file_name, header=None, parameters=None, data=None):
    """Save information as a pretty formatted text file.

    Parameters
    ----------
    file_name : str
        File name or path string to the desired save location. The ".txt"
        extension is not part of `file_name`.
    header : str
        Header string to be placed at the top of the file.
    parameters : List[List[str, str]]
        List of lists containing a parameter's label-value string pair.

        A label-value pair of `["Label String", "Value String"]` will be saved
        as "Label String: Value String" after the header.
    data : List[List[str, Quantity]]
        List of lists containing a data's label-value pair. The label is a
        string type and the value must be a `Quantity` object with a NumPy
        array `data` property.

        A label-value pair of `["Data Label", Quantity(np.array([1, 2]), u.m)]`
        will have the heading "Data Label (m)".
    """

    processed_file_name = file_name + ".txt"
    with open(processed_file_name, "w") as file:

        if header is not None:
            file.write(header)

        if (header is not None) and (parameters is not None or data is not None):
            file.write(BLANK_LINE_STRING)

        if parameters is not None:
            for parameter_count, [label, value] in enumerate(parameters):
                file.write(get_parameter_string(parameter_count,
                                                len(parameters), label, value))

        if (header is not None or parameters is not None) and (data is not None):
            file.write(BLANK_LINE_STRING)

        if data is not None:
            column_label_lengths = []

            column_header_string = ""
            column_underline_string = ""

            for data_count, [data_label, data_value] in enumerate(data):

                processed_column_label = data_label + " (" + \
                                         str(data_value.unit) + ")"
                column_label_lengths.append(len(processed_column_label))

                column_header_string += processed_column_label
                column_underline_string += "-" * len(processed_column_label)

                if data_count < len(data) - 1:
                    column_header_string += COLUMN_DELIMITER
                    column_underline_string += COLUMN_DELIMITER
                else:
                    column_header_string += NEWLINE_STRING
                    column_underline_string += NEWLINE_STRING

                reshaped_data_array = data_value.data.reshape((-1, 1))
                if data_count == 0:
                    raw_data = reshaped_data_array
                else:
                    raw_data = np.concatenate((raw_data, reshaped_data_array), axis=1)

            file.write(column_header_string)
            file.write(column_underline_string)

            number_of_rows, number_of_columns = raw_data.shape
            for row_index in range(number_of_rows):
                data_row_string = ""
                for column_index in range(number_of_columns):
                    data_value_string = str(raw_data[row_index][column_index])
                    data_row_string += data_value_string.ljust(column_label_lengths[column_index])

                    if column_index < number_of_columns - 1:
                        data_row_string += COLUMN_DELIMITER
                if row_index < number_of_rows - 1:
                    data_row_string += NEWLINE_STRING

                file.write(data_row_string)


def get_parameter_string(parameter_count, number_of_parameters, label, value):
    return label + ": " + value + (NEWLINE_STRING if parameter_count < number_of_parameters - 1 else "")
