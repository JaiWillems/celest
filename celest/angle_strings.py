

from celest.units.quantity import Quantity
from celest import units as u


DEGREE_SYMBOL = "\u00B0"
MINUTE_SYMBOL = "\u2032"
SECOND_SYMBOL = "\u2033"


def _ISO6709_representation(latitude: Quantity, longitude: Quantity, height:
                            Quantity) -> str:
    """Format geographical data in the ISO6709 standard.

    Parameters
    ----------
    latitude : Quantity
    longitude : Quantity
    height : Quantity

    Returns
    -------
    str
        ISO6709 standardized position string.
    """

    return " ".join([
        _get_latitude_string(latitude),
        _get_longitude_string(longitude),
        _get_height_string(height)
    ])


def _get_latitude_string(latitude: Quantity) -> str:
    sexigesimal_latitude = _get_sexagesimal_string(abs(latitude.to(u.deg)))
    direction = "N" if latitude.data >= 0 else "S"

    return sexigesimal_latitude + direction


def _get_sexagesimal_string(degree_angle: float) -> str:
    """Return the sexagesimal representation of a degree angle.

    It is assumed that the input angle is always positive.

    Parameters
    ----------
    degree_angle : float
        An absolute value degree angle.

    Returns
    -------
    str
        Sexagesimal representation string of `degree_angle`.
    """

    degree_string = "%.3s" % str(int(degree_angle)).zfill(2)
    minute_string = "%.2s" % str(int(degree_angle * 60 % 60)).zfill(2)
    second_string = str("%.2f" % (degree_angle * 3600 % 60)).zfill(5)

    return degree_string + DEGREE_SYMBOL + minute_string + MINUTE_SYMBOL + \
        second_string + SECOND_SYMBOL


def _get_longitude_string(longitude: Quantity) -> str:
    sexigesimal_longitude = _get_sexagesimal_string(abs(longitude.to(u.deg)))
    direction = "E" if longitude.data >= 0 else "W"

    return sexigesimal_longitude + direction


def _get_height_string(altitude: Quantity) -> str:
    return ("%.2f" % altitude.data) + str(altitude.unit)
