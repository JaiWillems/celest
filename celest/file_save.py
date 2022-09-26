

from celest.units.quantity import Quantity
from typing import IO, Any, Tuple
import numpy as np


NEWLINE_STRING = "\n"
BLANK_LINE_STRING = "\n\n"
COLUMN_DELIMITER = "\t"

HEADER_UNDERLINE_CHAR = "="
COLUMN_UNDERLINE_CHAR = "-"


class TextFileWriter:
    """TextFileWriter(file_name, header)

    Build and save data to a text file.

    Parameters
    ----------
    file_name : str
        File name or path to the desired save location.
    header : str
        Header to be placed at the top of the file.

    Methods
    -------
    add_layer(subheading, parameters, data, data_format="%.15f")
        Add a data layer to the text file.
    save()
        Save the text file.

    Examples
    --------
    Initialize a `TextFileWriter` object:

    >>> writer = TextFileWriter("file_name.txt", "Header String")

    Add a data layer to the text file:

    >>> writer.add_layer(
    ...     "Subheading String",
    ...     [["Label String", "Value String"]],
    ...     [["Data Label", Quantity(np.array([1, 2]), u.m)]],
    ...     "%.15f"
    ... )

    Save the data layers to a text file:

    >>> writer.save()
    """

    def __init__(self, file_name: str, header: str=None) -> None:
        self._file_name = _process_file_name(file_name)
        self._header = header
        self._layers = []

    def add_layer(self, subheading: str=None, parameters: list=None,
                  data: list=None, data_format: str="%.15f") -> None:
        """Add a data layer to the text file.

        A layer is a collection of information to be saved to a text file
        including a subheading, parameters, and data.

        Parameters
        ----------
        subheading : str
            Subheading to be placed at the top of the layer.
        parameters : list
            List of lists containing a parameter's label-value string pair.

            Each label-value pair of `["Label String", "Value String"]` will be
            saved as "Label String: Value String" after the subheading.
        data : list
            List of lists containing a data's label-value pair. The label is a
            string type and the value must be a `Quantity` object with a NumPy
            array `data` property or purely a NumPy array object.

            If the value is a `Quantity` object, the column header will be the
            name of the data label + " (" + the unit of the data + ")".
            Otherwise, it will just be the data label.
        data_format : str
            Data value format string.

        Examples
        --------
        Create the layer data:

        >>> subheading = "Subheading String"
        >>> parameters = [
        ...     ["Label 1", "Value 1"],
        ...     ["Label 2", "Value 2"]
        ... ]
        >>> data = [
        ...     ["Data Label", Quantity(np.array([1, 2]), u.m)],
        ...     ["Data Label 2", np.array([1, 2])]
        ... ]
        >>> data_format = "%.15f"

        Add the layer data to the text file:

        >>> writer = TextFileWriter("file_name.txt", "Header String")
        >>> writer.add_layer(subheading, parameters, data, data_format)
        """

        self._layers.append([subheading, parameters, data, data_format])

    def save(self) -> None:
        """Save stored data to a text file.

        This method creates a text file under the input file name and saves the
        header and layers to the file.
        """

        with open(self._file_name, "w", encoding='utf-8') as file:
            if self._header is not None:
                _write_header_to_file(file, self._header)

                if len(self._layers) > 0:
                    file.write(BLANK_LINE_STRING)

            for layer_count, layer_info in enumerate(self._layers):
                _write_layer_to_file(file, *layer_info)

                if layer_count < len(self._layers) - 1:
                    file.write(BLANK_LINE_STRING)


def _process_file_name(file_name: str) -> str:
    return file_name if file_name.endswith(".txt") else file_name + ".txt"


def _write_header_to_file(file: IO, header: str) -> None:
    file.write(header + NEWLINE_STRING + HEADER_UNDERLINE_CHAR * len(header))


def _write_layer_to_file(file: IO, subheading: str, parameters: list,
                         data: list, data_format: str) -> None:
    if subheading is not None:
        _write_subheading_to_file(file, subheading)
    if (subheading is not None) and (parameters is not None):
        file.write(BLANK_LINE_STRING)
    if parameters is not None:
        _write_parameters_to_file(file, parameters)
    if (subheading is not None or parameters is not None) and (data is not None):
        file.write(BLANK_LINE_STRING)
    if data is not None:
        _write_data_to_file(file, data, data_format)


def _write_subheading_to_file(file: IO, subheading: str) -> None:
    file.write(subheading)


def _write_parameters_to_file(file: IO, parameters: list) -> None:
    for parameter_count, [label, value] in enumerate(parameters):
        parameter_string = label + ": " + value
        if parameter_count < len(parameters) - 1:
            parameter_string += NEWLINE_STRING
        file.write(parameter_string)


def _write_data_to_file(file: IO, data: list, data_format: str) -> None:
    column_labels_string = ""
    column_underline_string = ""
    column_widths = []

    for data_count, [data_label, data_value] in enumerate(data):

        column_label, column_underline, column_width = _get_column_information(
            data_label,
            data_value,
            data_format
        )

        column_labels_string += column_label
        column_underline_string += column_underline
        column_widths.append(column_width)

        ending = COLUMN_DELIMITER if data_count < len(data) - 1 else NEWLINE_STRING
        column_labels_string += ending
        column_underline_string += ending

        column_data = _get_and_reshape_column_data(data_value)
        if data_count == 0:
            raw_data = column_data
        else:
            raw_data = np.concatenate((raw_data, column_data), axis=1)

    file.write(column_labels_string)
    file.write(column_underline_string)
    _write_column_data_to_file(file, column_widths, raw_data, data_format)


def _get_column_information(label: str, value: Any, value_format:
                            str) -> Tuple[str, str, int]:
    column_label = label
    if isinstance(value, Quantity):
        column_label += " (" + str(value.unit) + ")"

    value_width = _get_maximum_value_width(value, value_format)
    column_width = len(column_label)
    actual_width = max(value_width, column_width)

    processed_column_label = column_label.ljust(actual_width)
    column_underline = COLUMN_UNDERLINE_CHAR * actual_width

    return processed_column_label, column_underline, actual_width


def _get_maximum_value_width(value: Any, value_format: str) -> int:
    if isinstance(value, Quantity):
        return len(value_format % np.max(value.data))
    elif isinstance(value, np.ndarray):
        if isinstance(value.dtype, (int, float)):
            return len(value_format % np.max(value))
        else:
            return len(max(value.astype(str), key=len))
    else:
        return NotImplemented


def _get_and_reshape_column_data(value: Any) -> np.ndarray:
    if isinstance(value, Quantity):
        return value.data.reshape((-1, 1))
    else:
        return value.reshape((-1, 1))


def _write_column_data_to_file(file: IO, column_widths: list,
                               raw_data: np.ndarray, data_format: str) -> None:
    for row_count, row in enumerate(raw_data):
        row_string = ""
        for column_count, column in enumerate(row):
            row_string += _get_data_entry_string(
                column,
                column_widths[column_count],
                data_format
            )
            if column_count < len(row) - 1:
                row_string += COLUMN_DELIMITER
        if row_count < len(raw_data) - 1:
            row_string += NEWLINE_STRING
        file.write(row_string)


def _get_data_entry_string(value: Any, width: int, data_format: str) -> str:
    if isinstance(value, (int, float)):
        return (data_format % value).ljust(width)
    else:
        return str(value).ljust(width)
