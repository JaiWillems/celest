

import numpy as np


DEGREE_SYMBOL = "\u00B0"
MINUTE_SYMBOL = "\u2032"
SECOND_SYMBOL = "\u2033"


def sexagesimal(angles: np.ndarray) -> np.ndarray:
    """Convert decimal angles into sexagesimal angles.

    Parameters
    ----------
    angles : np.ndarray
        1-D array containing angles in decimal degrees.

    Returns
    -------
    np.ndarray
        1-D array containing sexagesimal angles as strings.

    Examples
    --------
    >>> angles = np.array([43.6532, -79.3832, -33.2833, 149.1000])
    >>> coor = Coordinate(...)
    >>> coor.sexagesimal(angles=angles)
    np.array(['+43°39′11.52″',
              '-79°22′59.52″',
              '-33°16′59.88″',
              '+149°06′00.00″'])
    """

    sexagesimal_angles = []
    for angle in angles:
        angle_sign = "+" if angle >= 0 else "-"
        sexagesimal_string = _get_sexagesimal_string(abs(angle))
        sexagesimal_angles.append(angle_sign + sexagesimal_string)

    return np.array(sexagesimal_angles)


def _get_sexagesimal_string(absolute_angle: float) -> str:

    degree_string = "%.3s" % str(int(absolute_angle)).zfill(2)
    minute_string = "%.2s" % str(int(absolute_angle * 60 % 60)).zfill(2)
    second_string = str("%.2f" % (absolute_angle * 3600 % 60)).zfill(5)

    return degree_string + DEGREE_SYMBOL + minute_string + MINUTE_SYMBOL + \
        second_string + SECOND_SYMBOL


def _ISO6709_representation(latitude: np.ndarray, longitude: np.ndarray,
                            altitude: np.ndarray) -> np.ndarray:
    """Format geographical data in the ISO6709 standard.

    Parameters
    ----------
    latitude : np.ndarray
        1-D array containing latitude angles in decimal degrees.
    longitude : np.ndarray
        1-D array containing longitude angles in decimal degrees.
    altitude : np.ndarray
        1-D array containing altitude angles in km.

    Returns
    -------
    np.ndarray
        1-D array containing ISO6709 standardized position strings.
    """

    iso_strings = []
    for lat, lon, alt in zip(latitude, longitude, altitude):
        latitude_string = _get_latitude_string(lat)
        longitude_string = _get_longitude_string(lon)
        altitude_string = _get_altitude_string(alt)

        iso_strings.append(latitude_string + " " + longitude_string + " " + altitude_string)

    return np.array(iso_strings)


def _get_latitude_string(latitude: float) -> str:

    sexigesimal_latitude = _get_sexagesimal_string(abs(latitude))
    direction = "N" if latitude >= 0 else "S"

    return sexigesimal_latitude + direction


def _get_longitude_string(longitude: float) -> str:

    sexigesimal_longitude = _get_sexagesimal_string(abs(longitude))
    direction = "E" if longitude >= 0 else "W"

    return sexigesimal_longitude + direction


def _get_altitude_string(altitude: float) -> str:

    return ("%.2f" % altitude) + "km"
