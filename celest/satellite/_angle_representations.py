"""Define various angle representations."""


from celest.core.decorators import set_module
import numpy as np


@set_module('celest.satellite')
def sexagesimal(angles: np.ndarray) -> np.ndarray:
    """Convert decimal angles into sexagesimal angles.

    Parameters
    ----------
    angles : np.ndarray
        Array of shape (n,) containing angles in decimal degrees.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing sexagesimal angles as strings.

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

    deg = u"\u00B0"
    min = u"\u2032"
    sec = u"\u2033"

    length = angles.shape[0]
    out_arr = np.empty((length,), dtype="<U32")

    for i in range(length):

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
    """Format geographical data in an international standard.

    This method formats the input geographical point location data in the
    ISO6709 standard format.

    Parameters
    ----------
    position : np.ndarray
        Array of shape (n, 3) with columns of geodetic latitude,
        terrestrial longitude, and geodetic altitude data.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing position strings in accordance to
        ISO6709 standard.
    """

    deg = "\u00B0"
    min = "\u2032"
    sec = "\u2033"

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
