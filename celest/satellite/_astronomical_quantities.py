

from statistics import mean
from typing import Tuple
import numpy as np


def nutation_angles(julian: np.ndarray) -> Tuple:
    """Return five Earth nutation angles.

    This method calculates the five nutation angles for each time in the
    `julian` array. The calculated angles include the mean elongation of
    the Moon from the Sun (D), mean anomaly of the Sun (M), mean anomaly of
    the Moon (N), Moon's argument of latitude (F), and the longitude of the
    ascending node of the Moon's mean orbit on the ecliptic measured from the
    mean equinox date (O).

    Parameters
    ----------
    jullian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    tuple
        Tuple of the form, `(D, M, N, F, O)`, where each item are 1-D NumPy
        arrays containing the Earth nutation angles in decimal degrees.

    Notes
    -----
    The time in Julian centeries since J2000, :math:`T`, can be calculated
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
    >>> julian = np.array([2446895.5])
    >>> D, M, N, F, O = nutation_angles(julian=julian)
    """

    T1, T2, T3 = _calculate_elapsed_JD_century_powers(julian, 3)

    D = _get_mean_elongation_moon_from_sun_deg(T1, T2, T3)
    M = _get_mean_anomaly_of_sun_deg(T1, T2, T3)
    N = _get_mean_anomaly_of_moon_deg(T1, T2, T3)
    F = _get_moons_argument_of_latitude_deg(T1, T2, T3)
    O = _get_longitude_of_ascending_node_deg(T1, T2, T3)

    return D, M, N, F, O


def _calculate_elapsed_JD_century_powers(julian, highest_power):

    elapsd_centuries = (julian - 2451545) / 36525
    return [elapsd_centuries ** i for i in range(1, highest_power + 1)]


def _get_mean_elongation_moon_from_sun_deg(T, T2, T3):

    return (297.85036 + 445267.111480 * T - 0.0019142 * T2 + T3 / 189474) % 360


def _get_mean_anomaly_of_sun_deg(T, T2, T3):

    return (357.52772 + 35999.050340 * T - 0.0001603 * T2 - T3 / 300000) % 360


def _get_mean_anomaly_of_moon_deg(T, T2, T3):

    return (134.96298 + 477198.867398 * T + 0.0086972 * T2 + T3 / 56250) % 360


def _get_moons_argument_of_latitude_deg(T, T2, T3):

    return (93.27191 + 483202.017538 * T - 0.0036825 * T2 + T3 / 327270) % 360


def _get_longitude_of_ascending_node_deg(T, T2, T3):

    return (125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000) % 360


def nutation_components(julian: np.ndarray) -> Tuple:
    """Return the nutations in longitude and obliquity.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    tuple
        Tuple of the form, `(longitude, obliquity)`, containing 1-D NumPy
        arrays with the nutation components of longitude and obliquity in
        decimal seconds.

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
    >>> julian = np.array([2446895.5])
    >>> delta_psi, delta_epsilon = nutation_componenets(julian=julian)
    """

    T1, T2, T3 = _calculate_elapsed_JD_century_powers(julian, 3)

    ascending_node_longitude = _get_longitude_of_ascending_node_deg(T1, T2, T3)
    sun_mean_lontitude = _get_mean_longitude_of_sun_deg(T1)
    moon_mean_longitude = _get_mean_longitude_of_moon_deg(T1)

    ascending_node_longitude = np.radians(ascending_node_longitude)
    sun_mean_lontitude = np.radians(sun_mean_lontitude)
    moon_mean_longitude = np.radians(moon_mean_longitude)

    nutation_in_longitude = _get_nutation_in_longitude(ascending_node_longitude,
                                                       sun_mean_lontitude,
                                                       moon_mean_longitude)
    nutation_in_obliquity = _get_nutation_in_obliquity(ascending_node_longitude,
                                                       sun_mean_lontitude,
                                                       moon_mean_longitude)

    return nutation_in_longitude, nutation_in_obliquity


def _get_mean_longitude_of_sun_deg(T):

    return (280.46645 + 36000.76983 * T)


def _get_mean_longitude_of_moon_deg(T):

    return 218.3165 + 481267.8813 * T


def _get_nutation_in_longitude(ascending_node_longitude, sun_mean_longitude,
                               moon_mean_longitude):

    P1 = -17.20 * np.sin(ascending_node_longitude)
    P2 = -1.32 * np.sin(2 * sun_mean_longitude)
    P3 = -0.23 * np.sin(2 * moon_mean_longitude)
    P4 = 0.21 * np.sin(2 * ascending_node_longitude)

    return P1 + P2 + P3 + P4


def _get_nutation_in_obliquity(ascending_node_longitude, sun_mean_longitude,
                               moon_mean_longitude):

    E1 = 9.20 * np.cos(ascending_node_longitude)
    E2 = 0.57 * np.cos(2 * sun_mean_longitude)
    E3 = 0.10 * np.cos(2 * moon_mean_longitude)
    E4 = -0.09 * np.cos(2 * ascending_node_longitude)

    return E1 + E2 + E3 + E4


def conventional_precession_angles(julian: np.ndarray) -> Tuple:
    """Return conventional precession angles.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    Tuple
        Tuple containing precession angles in decimal arcseconds.

    Notes
    -----
    The convetional precession angles are those derived using the IAU 2000A
    model. [SL13c]_

    References
    ----------
    .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
       Astronomy and Astrophysics Library. Springer-Verlag, 2013, pp. 219.
    """

    T1, T2, T3, T4, T5 = _calculate_elapsed_JD_century_powers(julian, 5)

    zeta = _calculate_zeta_precession_angle_arcseconds(T1, T2, T3, T4, T5)
    theta = _calculate_theta_precession_angle_arcseconds(T1, T2, T3, T4, T5)
    z = _calculate_z_precession_angle_arcseconds(T1, T2, T3, T4, T5)

    return zeta, theta, z


def _calculate_zeta_precession_angle_arcseconds(T1, T2, T3, T4, T5):

    return 2.59796176 + 2306.0809506 * T1 + 0.3019015 * T2 + 0.0179663 * T3 - \
        0.0000327 * T4 - 0.0000002 * T5


def _calculate_theta_precession_angle_arcseconds(T1, T2, T3, T4, T5):

    return 2004.1917476 * T1 - 0.4269353 * T2 - 0.0418251 * T3 - \
        0.0000601 * T4 - 0.0000001 * T5


def _calculate_z_precession_angle_arcseconds(T1, T2, T3, T4, T5):

    return - 2.5976176 + 2306.0803226 * T1 + 1.0947790 * T2 + 0.0182273 * T3 + \
        0.0000470 * T4 - 0.0000003 * T5


def mean_obliquity(julian: np.ndarray) -> np.ndarray:
    """Return the ecliptic's mean obliquity in the J2000 epoch.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        1-D array containing the ecliptic's mean obliquity in degrees and
        decimals.

    See Also
    --------
    apparent_obliquity :
        Return the ecliptic's apparent obliquity in the J2000 epoch.

    Notes
    -----
    The methods to calculate the mean obliquity of the ecliptic are given in
    "Astronomical Algorithms" by Jean Meeus and are only valid for 10000 years
    on either side of J2000. [Mee98d]_

    References
    ----------
    .. [Mee98d] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 147. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> epsilon_0 = mean_obliquity(julian=julian)
    """

    T = _calculate_elapsed_JD_century_powers(julian, 1)[0]
    U = T / 100

    mean_obliquity = (84381.448 - 4680.93 * U - 1.55 * U ** 2 +
                      1999.25 * U ** 3 - 51.38 * U ** 4 - 249.67 * U ** 5 -
                      39.05 * U ** 6 + 7.12 * U ** 7 + 27.87 * U ** 8 +
                      5.79 * U ** 9 + 2.45 * U ** 10) / 3600

    return mean_obliquity


def apparent_obliquity(julian: np.ndarray) -> np.ndarray:
    """Return the ecliptic's apparent obliquity in the J2000 epoch.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        1-D array containing the ecliptic's apparent obliquity in degrees and
        decimals.

    See Also
    --------
    mean_obliquity :
        Return the ecliptic's mean obliquity in the J2000 epoch.

    Notes
    -----
    The methods to calculate the apparent obliquity of the ecliptic are given
    in "Astronomical Algorithms" by Jean Meeus and are only valid for 10000
    years on either side of J2000. [Mee98e]_

    References
    ----------
    .. [Mee98e] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 147. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> epsilon = apparent_obliquity(julian=julian)
    """

    average_obliquity = mean_obliquity(julian)
    _, nutation_in_obliquity = nutation_components(julian)

    return average_obliquity + nutation_in_obliquity / 3600


def from_julian(julian: np.ndarray) -> Tuple:
    """Return year, month, day from Julian times.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    Tuple
        Returns a tuple, `(year, month, day)`, containing 1-D NumPy arrays
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
    >>> julian = np.array([2436116.31])
    >>> apparent_obliquity(julian=julian)
    (1957, 10, 4.81)
    """

    JD = julian + 0.5
    A, F = JD.astype(int), JD % 1

    ind = np.where(A >= 2291161)[0]
    alpha = ((A[ind] - 1867216.25) / 36524.25).astype(int)
    A[ind] = A + 1 + alpha - (alpha / 4).astype(int)

    B = A + 1524
    C = ((B - 122.1) / 365.25).astype(int)
    D = (365.25 * C).astype(int)
    E = ((B - D) / 30.6001).astype(int)

    day = B - D - (30.6001 * E).astype(int) + F

    month = E
    idx1 = np.where(month < 14)[0]
    idx2 = np.where((month == 14) | (month == 15))[0]
    month[idx1] = month[idx1] - 1
    month[idx2] = month[idx2] - 13

    year = C
    idx1 = np.where(month > 2)[0]
    idx2 = np.where((month == 1) | (month == 2))[0]
    year[idx1] = year[idx1] - 4716
    year[idx2] = year[idx2] - 4715

    return year, month, day


def day_of_year(julian: np.ndarray) -> np.ndarray:
    """Return the day of the year.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

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
    >>> julData = np.array([2447273.5])
    >>> AstronomicalQuantities().day_of_year(julData=julData)
    113
    """

    year, month, day = from_julian(julian)

    K = np.full(day.shape, 2)
    K[np.where(year % 4 == 0)] = 1

    N = (275 * month / 9).astype(int)
    N = N - K * ((month + 9) / 12).astype(int) + day - 30

    return N.astype(int)


def equation_of_time(julian: np.ndarray) -> np.ndarray:
    """Return the equation of time in decimal degrees.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        1-D array containing the Equation of Time in decimal degrees.

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
    >>> julian = np.array([2448908.5])
    >>> equation_of_time(julian=julian)
    3.427351
    """

    T1, T2, T3 = _calculate_elapsed_JD_century_powers(julian, 3)

    sun_mean_longitude = np.radians(_get_mean_longitude_of_sun_deg(T1))
    sun_mean_anomaly = np.radians(_get_mean_anomaly_of_sun_deg(T1, T2, T3))
    earth_eccentricity = _get_earth_eccentricity(T1)

    sun_apparent_obliquity = np.radians(apparent_obliquity(julian))
    y = np.tan(sun_apparent_obliquity / 2) ** 2

    equation_of_time = y * np.sin(2 * sun_mean_longitude) - \
        2 * earth_eccentricity * np.sin(sun_mean_anomaly) + \
        4 * earth_eccentricity * y * np.sin(sun_mean_anomaly) * np.cos(2 * sun_mean_longitude) - \
        0.5 * y ** 2 * np.sin(4 * sun_mean_longitude) - \
        1.25 * earth_eccentricity ** 2 * np.sin(2 * sun_mean_anomaly)

    return np.degrees(equation_of_time)


def _get_earth_eccentricity(T):

    return 0.016708634 - 0.000042037 * T - 0.0000001267 * T ** 2


def equation_of_equinoxes(julian: np.ndarray) -> np.ndarray:
    """Return the equation of the equinoxes in decimal arcseconds.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        1-D array containing the equation of equinoxes in decimal arcseconds.

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
    >>> julian = np.array([2446895.5]])
    >>> equation_of_equinoxes(julian=julian)
    -0.2317
    """

    nutation_in_longitude, _ = nutation_components(julian)
    obliquity = np.radians(apparent_obliquity(julian))

    return nutation_in_longitude * np.cos(obliquity) / 15


def sun_right_ascension(julian: np.ndarray) -> np.ndarray:
    """Return the right ascension of the mean sun position.

    This method calculates the Sun's right ascension to an accuracy of 0.01 of
    a degree.

    Parameters
    ----------
    julian : np.ndarray
        1-D array containing Julian times in the J2000 epoch.

    Returns
    -------
    np.ndarray
        1-D array containing containing the right ascension of the mean sun in
        decimal degrees.

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
    >>> julian = np.array([2448908.5])
    >>> sun_right_ascension(julian=julian)
    np.array([198.38166])
    """

    T1, T2, T3 = _calculate_elapsed_JD_century_powers(julian, 3)

    mean_sun_longitude = _get_mean_longitude_of_sun_deg(T1)
    mean_sun_anomaly = np.radians(_get_mean_anomaly_of_sun_deg(T1, T2, T3))

    sun_center = (1.914602 - 0.004817 * T1 - 0.000014 * T2) * np.sin(mean_sun_anomaly) + \
        (0.019993 - 0.000101 * T1) * np.sin(2 * mean_sun_anomaly) + \
        0.000289 * np.sin(3 * mean_sun_anomaly)

    sun_true_longitude = np.radians(mean_sun_longitude + sun_center)
    obliquity = np.radians(mean_obliquity(julian))

    alpha = np.arctan2(np.cos(obliquity) * np.sin(sun_true_longitude), np.cos(sun_true_longitude))

    return np.degrees(alpha) % 360
