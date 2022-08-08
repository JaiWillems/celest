

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
        Latitude angle in decimal degrees.
    longitude : Quantity
        Longitude angle in decimal degrees.
    height : Quantity
        Altitude in km.

    Returns
    -------
    str
        ISO6709 standardized position string.
    """

    latitude_string = _get_latitude_string(latitude)
    longitude_string = _get_longitude_string(longitude)
    altitude_string = _get_altitude_string(height)

    return latitude_string + " " + longitude_string + " " + altitude_string


def _get_latitude_string(latitude: Quantity) -> str:
    sexigesimal_latitude = _get_sexagesimal_string(abs(latitude.to(u.deg)))
    direction = "N" if latitude.data >= 0 else "S"

    return sexigesimal_latitude + direction


def _get_sexagesimal_string(absolute_deg_angle: float) -> str:
    degree_string = "%.3s" % str(int(absolute_deg_angle)).zfill(2)
    minute_string = "%.2s" % str(int(absolute_deg_angle * 60 % 60)).zfill(2)
    second_string = str("%.2f" % (absolute_deg_angle * 3600 % 60)).zfill(5)

    return degree_string + DEGREE_SYMBOL + minute_string + MINUTE_SYMBOL + \
        second_string + SECOND_SYMBOL


def _get_longitude_string(longitude: Quantity) -> str:
    sexigesimal_longitude = _get_sexagesimal_string(abs(longitude.to(u.deg)))
    direction = "E" if longitude.data >= 0 else "W"

    return sexigesimal_longitude + direction


def _get_altitude_string(altitude: Quantity) -> str:
    return ("%.2f" % altitude.data) + str(altitude.unit)
