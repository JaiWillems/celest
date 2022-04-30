

import numpy as np


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

    deg, min, sec = u"\u00B0", u"\u2032", u"\u2033"

    n = angles.shape[0]
    out_arr = np.empty((n,), dtype="<U32")

    for i in range(n):

        ang = angles[i]

        sign = "+" if ang >= 0 else "-"
        degree = int(abs(ang))
        minute = int(60 * np.round(abs(ang) - degree, 2))
        second = np.round(abs(60 * (60 * (abs(ang) - degree) - minute)), 2)

        degree = "%.3s" % str(degree).zfill(2)
        minute = "%.2s" % str(minute).zfill(2)
        second = str("%.2f" % second).zfill(5)

        out_arr[i] = f"{sign}{degree}{deg}{minute}{min}{second}{sec}"

    return out_arr


def _ISO6709_representation(position: np.ndarray) -> np.ndarray:
    """Format geographical data in the ISO6709 standard.

    Parameters
    ----------
    position : np.ndarray
        2-D array containing columns of latitude, longitude, and altitude data.

    Returns
    -------
    np.ndarray
        1-D array containing ISO6709 standardized position strings.
    """

    deg, min, sec = "\u00B0", "\u2032", "\u2033"

    out_arr = np.empty((position.shape[0],), dtype="<U37")

    for i in range(position.shape[0]):

        lat, lon, h = position[i, 0], position[i, 1], position[i, 2]

        # Format latitude string.
        degree = int(abs(lat))
        minute = int(60 * np.round(abs(lat) - degree, 2))
        second = np.round(abs(60 * (60 * (abs(lat) - degree) - minute)), 2)
        direction = "N" if lat >= 0 else "S"

        degree = "%.2s" % str(degree).zfill(2)
        minute = "%.2s" % str(minute).zfill(2)
        second = str("%.2f" % second).zfill(5)

        lat_str = f"{degree}{deg}{minute}{min}{second}{sec}{direction} "

        # Format longitude string.
        degree = int(abs(lon))
        minute = int(60 * np.round(abs(lon) - degree, 2))
        second = np.round(abs(60 * (60 * (abs(lon) - degree) - minute)), 2)
        direction = "E" if lon >= 0 else "W"

        degree = "%.3s" % str(degree).zfill(2)
        minute = "%.2s" % str(minute).zfill(2)
        second = str("%.2f" % second).zfill(5)

        lon_str = f"{degree}{deg}{minute}{min}{second}{sec}{direction} "

        # Format height string.
        h = "%.2f" % h
        h_str = f"{h}km"

        out_arr[i] = lat_str + lon_str + h_str

    return out_arr
