

from celest.constants import (
    JD2000_DATE,
    DAYS_IN_JULIAN_CENTURY,
    DAYS_PER_YEAR,
    EARTH_ROTATION_RATE_DEG_PER_DAY,
    ERA_AT_JD2000_DEG,
    WGS84_MAJOR_AXIS_KM,
    WGS83_FIRST_ECCENTRICITY
)
from celest.units.quantity import Quantity
from celest import units as u
from math import sin, sqrt
from typing import Tuple
import numpy as np


def earth_rotation_angle(julian: Quantity) -> Quantity:
    """Return Earth rotation angle.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Quantity
        Earth rotation angle.

    Notes
    -----
    The Earth rotation angle is calculated using the following:

    .. math:: \gamma^\circ = 360.9856123035484\Delta T + 280.46

    where :math:`\Delta T=JD-2451545` is the elapsed days since the JD2000
    epoch where :math:`JD` is the Julian day. [Kok17a]_

    References
    ----------
    .. [Kok17a] Don Koks. Changing Coordinates in the Context of Orbital
       Mechanics. Cyber and Electronic Warfare Division, Defence Science,
       and Technology Group, Jan.2017, p. 12 - 13.

    Examples
    --------
    >>> julian = Quantity(30462.50000, u.jd2000)
    >>> era = earth_rotation_angle(julian)
    """

    days_since_jd2000 = julian.to(u.jd2000) - JD2000_DATE
    earth_rotation_angles = (EARTH_ROTATION_RATE_DEG_PER_DAY *
                             days_since_jd2000 + ERA_AT_JD2000_DEG) % 360
    return Quantity(earth_rotation_angles, u.deg)


def nutation_angles(julian: Quantity) -> Tuple:
    """Return five Earth nutation angles.

    This method calculates the five nutation angles for each time in the
    `julian` array. The calculated angles include the mean elongation of
    the Moon from the Sun (D), mean anomaly of the Sun (M), mean anomaly of
    the Moon (N), Moon's argument of latitude (F), and the longitude of the
    ascending node of the Moon's mean orbit on the ecliptic measured from the
    mean equinox date (O).

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    tuple
        Tuple of the form, `(D, M, N, F, O)`, where each item is a Quantity
        object containing the Earth nutation angles.

    Notes
    -----
    The time in Julian centeries since JD2000, :math:`T`, can be calculated
    from the Julian day, :math:`JD`, from the following:

    .. math:: T = \\frac{JD - 2451545}{36525}

    The nutation angles can then be calculated in decimal degrees. [Mee98b]_

    .. math:: D = 297.85036 + 445267.111480T - 0.0019142T^2 + T^3 / 189474

    .. math:: M = 357.52772 + 35999.050340T - 0.0001603T^2 - T^3 / 300000

    .. math:: N = 134.96298 + 477198.867398T + 0.0086972T^2 + T^3 / 56250

    .. math:: F = 93.27191 + 483202.017538T - 0.0036825T^2 + T^3 / 327270

    .. math:: O = 125.04452 - 1934.136261T + 0.0020708T^2 + T^3 / 450000

    References
    ----------
    .. [Mee98b] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 143-144. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2446895.5, u.jd2000)
    >>> D, M, N, F, O = nutation_angles(julian=julian)
    """

    return (
        get_mean_elongation_moon_from_sun_deg(julian),
        get_mean_anomaly_of_sun(julian),
        get_mean_anomaly_of_moon_deg(julian),
        get_moons_argument_of_latitude_deg(julian),
        get_longitude_of_ascending_node_deg(julian)
    )


def get_mean_elongation_moon_from_sun_deg(julian: Quantity) -> Quantity:
    t1, t2, t3 = _calculate_raw_elapsed_jd_century_powers(julian, 3)
    result = (297.85036 + 445267.111480 * t1 - 0.0019142 * t2 + t3 / 189474) % 360
    return Quantity(result, u.deg)


def _calculate_raw_elapsed_jd_century_powers(julian: Quantity, highest_power):
    centuries_since_jd2000 = (julian.to(u.jd2000) - JD2000_DATE) / DAYS_IN_JULIAN_CENTURY
    return [centuries_since_jd2000 ** p for p in range(1, highest_power + 1)]


def get_mean_anomaly_of_sun(julian: Quantity) -> Quantity:
    t1, t2, t3 = _calculate_raw_elapsed_jd_century_powers(julian, 3)
    result = (357.52772 + 35999.050340 * t1 - 0.0001603 * t2 - t3 / 300000) % 360
    return Quantity(result, u.deg)


def get_mean_anomaly_of_moon_deg(julian: Quantity) -> Quantity:
    t1, t2, t3 = _calculate_raw_elapsed_jd_century_powers(julian, 3)
    result = (134.96298 + 477198.867398 * t1 + 0.0086972 * t2 + t3 / 56250) % 360
    return Quantity(result, u.deg)


def get_moons_argument_of_latitude_deg(julian: Quantity) -> Quantity:
    t1, t2, t3 = _calculate_raw_elapsed_jd_century_powers(julian, 3)
    result = (93.27191 + 483202.017538 * t1 - 0.0036825 * t2 + t3 / 327270) % 360
    return Quantity(result, u.deg)


def get_longitude_of_ascending_node_deg(julian):
    t1, t2, t3 = _calculate_raw_elapsed_jd_century_powers(julian, 3)
    result = (125.04452 - 1934.136261 * t1 + 0.0020708 * t2 + t3 / 450000) % 360
    return Quantity(result, u.deg)


def nutation_components(julian: Quantity) -> Tuple:
    """Return the nutations in longitude and obliquity.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    tuple
        Tuple of the form, `(longitude, obliquity)`, containing Quantity
        objects with the nutation components of longitude and obliquity.

    Notes
    -----
    The methods to calculate the nutations in longitude and obliquity are taken
    from "Astronomical Algorithms" by Jean Meeus. [Mee98c]_

    References
    ----------
    .. [Mee98c] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 143-144. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2446895.5, u.jd2000)
    >>> delta_psi, delta_epsilon = nutation_componenets(julian=julian)
    """

    ascending_node_longitude = get_longitude_of_ascending_node_deg(julian)
    sun_mean_longitude = get_mean_longitude_of_sun(julian)
    moon_mean_longitude = get_mean_longitude_of_moon(julian)

    nutation_in_longitude = get_nutation_in_longitude(ascending_node_longitude,
                                                      sun_mean_longitude,
                                                      moon_mean_longitude)
    nutation_in_obliquity = get_nutation_in_obliquity(ascending_node_longitude,
                                                      sun_mean_longitude,
                                                      moon_mean_longitude)

    return nutation_in_longitude, nutation_in_obliquity


def get_mean_longitude_of_sun(julian: Quantity) -> Quantity:
    t, = _calculate_raw_elapsed_jd_century_powers(julian, 1)
    return Quantity(280.46645 + 36000.76983 * t, u.deg)


def get_mean_longitude_of_moon(julian: Quantity) -> Quantity:
    t, = _calculate_raw_elapsed_jd_century_powers(julian, 1)
    return Quantity(218.3165 + 481267.8813 * t, u.deg)


def get_nutation_in_longitude(ascending_node_longitude: Quantity,
                              sun_mean_longitude: Quantity,
                              moon_mean_longitude: Quantity) -> Quantity:
    p1 = -17.20 * np.sin(ascending_node_longitude.to(u.rad))
    p2 = -1.32 * np.sin(2 * sun_mean_longitude.to(u.rad))
    p3 = -0.23 * np.sin(2 * moon_mean_longitude.to(u.rad))
    p4 = 0.21 * np.sin(2 * ascending_node_longitude.to(u.rad))

    return Quantity(p1 + p2 + p3 + p4, u.arcsec)


def get_nutation_in_obliquity(ascending_node_longitude: Quantity,
                              sun_mean_longitude: Quantity,
                              moon_mean_longitude: Quantity) -> Quantity:
    e1 = 9.20 * np.cos(ascending_node_longitude.to(u.rad))
    e2 = 0.57 * np.cos(2 * sun_mean_longitude.to(u.rad))
    e3 = 0.10 * np.cos(2 * moon_mean_longitude.to(u.rad))
    e4 = -0.09 * np.cos(2 * ascending_node_longitude.to(u.rad))

    return Quantity(e1 + e2 + e3 + e4, u.arcsec)


def conventional_precession_angles(julian: Quantity) -> Tuple:
    """Return conventional precession angles.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Tuple
        Precession angles.

    Notes
    -----
    The conventional precession angles are those derived using the IAU 2000A
    model. [SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 219.
    """

    return (
        _get_zeta_precession_angle(julian),
        _get_theta_precession_angle(julian),
        _get_z_precession_angle(julian)
    )


def _get_zeta_precession_angle(julian: Quantity) -> Quantity:
    t1, t2, t3, t4, t5 = _calculate_raw_elapsed_jd_century_powers(julian, 5)
    result = 2.59796176 + 2306.0809506 * t1 + 0.3019015 * t2 + \
        0.0179663 * t3 - 0.0000327 * t4 - 0.0000002 * t5
    return Quantity(result, u.arcsec)


def _get_theta_precession_angle(julian: Quantity) -> Quantity:
    t1, t2, t3, t4, t5 = _calculate_raw_elapsed_jd_century_powers(julian, 5)
    result = 2004.1917476 * t1 - 0.4269353 * t2 - 0.0418251 * t3 - \
        0.0000601 * t4 - 0.0000001 * t5
    return Quantity(result, u.arcsec)


def _get_z_precession_angle(julian: Quantity) -> Quantity:
    t1, t2, t3, t4, t5 = _calculate_raw_elapsed_jd_century_powers(julian, 5)
    result = - 2.5976176 + 2306.0803226 * t1 + 1.0947790 * t2 + \
        0.0182273 * t3 + 0.0000470 * t4 - 0.0000003 * t5
    return Quantity(result, u.arcsec)


def mean_obliquity(julian: Quantity) -> Quantity:
    """Return the ecliptic's mean obliquity in the JD2000 epoch.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Quantity
        The ecliptic's mean obliquity.

    See Also
    --------
    apparent_obliquity :
        Return the ecliptic's apparent obliquity in the JD2000 epoch.

    Notes
    -----
    The methods to calculate the mean obliquity of the ecliptic are given in
    "Astronomical Algorithms" by Jean Meeus and are only valid for 10000 years
    on either side of JD2000. [Mee98d]_

    References
    ----------
    .. [Mee98d] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 147. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2446895.5, u.jd2000)
    >>> epsilon_0 = mean_obliquity(julian=julian)
    """

    t, = _calculate_raw_elapsed_jd_century_powers(julian, 1)
    U = t / 100

    result = (84381.448 - 4680.93 * U - 1.55 * U ** 2 + 1999.25 * U ** 3 -
              51.38 * U ** 4 - 249.67 * U ** 5 - 39.05 * U ** 6 + 7.12 * U **
              7 + 27.87 * U ** 8 + 5.79 * U ** 9 + 2.45 * U ** 10) / 3600

    return Quantity(result, u.deg)


def apparent_obliquity(julian: Quantity) -> Quantity:
    """Return the ecliptic's apparent obliquity in the JD2000 epoch.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Quantity
        The ecliptic's apparent obliquity.

    See Also
    --------
    mean_obliquity :
        Return the ecliptic's mean obliquity in the JD2000 epoch.

    Notes
    -----
    The methods to calculate the apparent obliquity of the ecliptic are given
    in "Astronomical Algorithms" by Jean Meeus and are only valid for 10000
    years on either side of JD2000. [Mee98e]_

    References
    ----------
    .. [Mee98e] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 147. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2446895.5, u.jd2000)
    >>> epsilon = apparent_obliquity(julian=julian)
    """

    average_obliquity = mean_obliquity(julian)
    _, nutation_in_obliquity = nutation_components(julian)

    return average_obliquity + nutation_in_obliquity


def from_julian(julian: Quantity) -> Tuple:
    """Return year, month, day from Julian times.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Tuple
        Tuple of the form, `(year, month, day)`, containing 1-D NumPy arrays
        representing the year, month, and decimal days of the input times.

    Notes
    -----
    Refer to page 63 of Astronomical Algorithms by Jean Meeus for a detailed
    implementation of the algorithm. [Mee98f]_

    References
    ----------
    .. [Mee98f] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 63. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2436116.31, u.jd2000)
    >>> apparent_obliquity(julian=julian)
    (1957, 10, 4.81)
    """

    if isinstance(julian.to(u.jd2000), (float, int)):
        jd = np.array([julian.to(u.jd2000)]) + 0.5
    else:
        jd = np.array(julian.to(u.jd2000)) + 0.5
    a, f = jd.astype(int), jd % 1

    ind = np.where(a >= 2291161)[0]
    alpha = ((a[ind] - 1867216.25) / DAYS_IN_JULIAN_CENTURY).astype(int)
    a[ind] = a + 1 + alpha - (alpha / 4).astype(int)

    b = a + 1524
    c = ((b - 122.1) / DAYS_PER_YEAR).astype(int)
    d = (DAYS_PER_YEAR * c).astype(int)
    e = ((b - d) / 30.6001).astype(int)

    day = b - d - (30.6001 * e).astype(int) + f

    month = e
    idx1 = np.where(month < 14)[0]
    idx2 = np.where((month == 14) | (month == 15))[0]
    month[idx1] = month[idx1] - 1
    month[idx2] = month[idx2] - 13

    year = c
    idx1 = np.where(month > 2)[0]
    idx2 = np.where((month == 1) | (month == 2))[0]
    year[idx1] = year[idx1] - 4716
    year[idx2] = year[idx2] - 4715

    return year, month, day


def day_of_year(julian: Quantity) -> np.ndarray:
    """Return the day of the year.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    np.ndarray
        1-D array containing the day of the year.

    Notes
    -----
    The day of year is calculated using the methods in "Astronomical
    Algorithms" by Jean Meeus. [Mee98g]_

    References
    ----------
    .. [Mee98g] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 65. isbn: 9780943396613.

    Examples
    --------
    >>> julData = Quantity(2447273.5, u.jd2000)
    >>> day_of_year(julian=julian)
    113
    """

    year, month, day = from_julian(julian)

    K = np.full(day.shape, 2)
    K[np.where(year % 4 == 0)] = 1

    N = (275 * month / 9).astype(int)
    N = N - K * ((month + 9) / 12).astype(int) + day - 30

    return N.astype(int)


def equation_of_time(julian: Quantity) -> Quantity:
    """Return the equation of time.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Quantity
        Equation of Time.

    Notes
    -----
    The equation of time is calculated using the methods described in
    "Astronomical Algorithms" by Jean Meeus. [Mee98h]_

    References
    ----------
    .. [Mee98h] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 184 - 185. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2448908.5, u.jd2000)
    >>> equation_of_time(julian=julian)
    3.427351
    """

    sun_mean_longitude = get_mean_longitude_of_sun(julian).to(u.rad)
    sun_mean_anomaly = get_mean_anomaly_of_sun(julian).to(u.rad)
    earth_eccentricity = get_earth_eccentricity(julian).to(u.rad)

    sun_apparent_obliquity = apparent_obliquity(julian).to(u.rad)
    y = np.tan(sun_apparent_obliquity / 2) ** 2

    e1 = y * np.sin(2 * sun_mean_longitude)
    e2 = - 2 * earth_eccentricity * np.sin(sun_mean_anomaly)
    e3 = 4 * earth_eccentricity * y * np.sin(sun_mean_anomaly) * np.cos(2 * sun_mean_longitude)
    e4 = - 0.5 * y ** 2 * np.sin(4 * sun_mean_longitude)
    e5 = - 1.25 * earth_eccentricity ** 2 * np.sin(2 * sun_mean_anomaly)

    return Quantity(e1 + e2 + e3 + e4 + e5, u.rad)


def get_earth_eccentricity(julian: Quantity) -> Quantity:
    t, = _calculate_raw_elapsed_jd_century_powers(julian, 1)
    result = 0.016708634 - 0.000042037 * t - 0.0000001267 * t ** 2
    return Quantity(result, u.rad)


def equation_of_equinoxes(julian: Quantity) -> Quantity:
    """Return the equation of the equinoxes.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Quantity
        Equation of equinoxes.

    Notes
    -----
    The equation of the equinoxes is given by
    :math:`\frac{1}{15}\Delta\psi\cos\epsilon` where :math:`\Delta\psi` is the
    nutation in longitude represented in arcseconds and, :math:`\epsilon` is
    the true obliquity of the ecliptic. [Mee98i]_

    References
    ----------
    .. [Mee98i] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 88. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2446895.5, u.jd2000)
    >>> equation_of_equinoxes(julian=julian)
    -0.2317
    """

    nutation_in_longitude, _ = nutation_components(julian)
    obliquity = apparent_obliquity(julian).to(u.rad)
    result = nutation_in_longitude.data * np.cos(obliquity) / 15

    return Quantity(result, u.arcsec)


def sun_right_ascension(julian: Quantity) -> Quantity:
    """Return the right ascension of the mean sun position.

    This method calculates the Sun's right ascension to an accuracy of 0.01 of
    a degree.

    Parameters
    ----------
    julian : Quantity
        Quantity object containing date data in the JD2000 epoch.

    Returns
    -------
    Quantity
        Right ascension of the mean sun position.

    Notes
    -----
    The right ascension of the Sun is calculated using the methods discussed in
    "Astronomical Algorithms" by Jean Meeus. [Mee98j]_

    References
    ----------
    .. [Mee98j] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 163 - 165. isbn: 9780943396613.

    Examples
    --------
    >>> julian = Quantity(2448908.5, u.jd2000)
    >>> sun_right_ascension(julian=julian)
    np.array([198.38166])
    """

    mean_sun_longitude = get_mean_longitude_of_sun(julian).to(u.deg)
    mean_sun_anomaly = get_mean_anomaly_of_sun(julian).to(u.rad)

    t1, t2 = _calculate_raw_elapsed_jd_century_powers(julian, 2)
    sun_center = (1.914602 - 0.004817 * t1 - 0.000014 * t2) * \
        np.sin(mean_sun_anomaly) + (0.019993 - 0.000101 * t1) * \
        np.sin(2 * mean_sun_anomaly) + 0.000289 * np.sin(3 * mean_sun_anomaly)

    sun_true_longitude = np.radians(mean_sun_longitude + sun_center)
    obliquity = mean_obliquity(julian).to(u.rad)

    alpha = np.arctan2(np.cos(obliquity) * np.sin(sun_true_longitude),
                       np.cos(sun_true_longitude))

    return Quantity(alpha % (2 * np.pi), u.rad)


def radius_of_curvature_in_prime_vertical(latitude: Quantity) -> Quantity:
    return Quantity(WGS84_MAJOR_AXIS_KM / sqrt(1 -
                    WGS83_FIRST_ECCENTRICITY ** 2 * sin(latitude.to(u.rad)) ** 2), u.km)
