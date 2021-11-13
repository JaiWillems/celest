"""Astronomical quantities for astronomical calculations.

This module contains the functionality to cary out advanced astronomical
calculations.
"""


from typing import Tuple
import numpy as np


def nutation_angles(julian: np.ndarray) -> Tuple:
    """Return the five Earth nutation angles.

    This method takes a Julian day array and calculates the five nutation
    angles at each time. The calculated angles include the mean elongation of
    the Moon from the Sun (D), mean anomaly of the Sun (M), mean anomaly of
    the Moon (N), Moon's argument of latitude (F), and the longitude of the
    ascending node of the Moon's mean orbit on the ecliptic measured from the
    mean equinox date (O).

    Parameters
    ----------
    jullian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    tuple
        Returns a tuple, `(D, M, N, F, O)`, where the items are Numpy arrays
        of shape (n,) containing the Earth nutation angles in degrees and
        decimals.

    Notes
    -----
    The time in Julian centeries since J2000, :math:`T`, can be calculated
    from the Julian day, :math:`JD`, from the following:

    .. math:: T = \\frac{JD - 2451545}{36525}

    The nutation angles can then be calculated in degrees and decimals. [Mee98b]_

    .. math:: D = 297.85036 + 445267.111480T - 0.0019142T^2 + T^3 / 189474

    .. math:: M = 357.52772 + 35999.050340T - 0.0001603T^2 - T^3 / 300000

    .. math:: N = 134.96298 + 477198.867398T + 0.0086972T^2 + T^3 / 56250

    .. math:: F = 93.27191 + 483202.017538T - 0.0036825T^2 + T^3 / 327270

    .. math:: O = 125.04452 - 1934.136261T + 0.0020708T^2 + T^3 / 450000

    References
    ----------
    .. [Mee98b] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 143–144. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> D, M, N, F, O = nutation_angles(julian=julian)
    """

    T = (julian - 2451545) / 36525
    T2, T3 = T ** 2, T ** 3

    D = (297.85036 + 445267.111480 * T - 0.0019142 * T2 + T3 / 189474) % 360
    M = (357.52772 + 35999.050340 * T - 0.0001603 * T2 - T3 / 300000) % 360
    N = (134.96298 + 477198.867398 * T + 0.0086972 * T2 + T3 / 56250) % 360
    F = (93.27191 + 483202.017538 * T - 0.0036825 * T2 + T3 / 327270) % 360
    O = (125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000) % 360

    return D, M, N, F, O

def nutation_components(julian: np.ndarray) -> Tuple:
    """Return the nutations in longitude and obliquity.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    tuple
        Returns a tuple, `(longitude, obliquity)`, where the items are NumPy
        arrays of shape (n,) containing the nutation components of longitude
        and obliquity in seconds and decimals.

    Notes
    -----
    The time in Julian centeries since J2000, :math:`T`, can be calculated
    from the Julian day, :math:`JD`, from the following:

    .. math:: T = \\frac{JD - 2451545}{36525}

    The longitude of the ascending node of the Moon's orbit on the ecliptic
    measured from the mean equinox can be calculated in decimal degrees
    using the following:

    .. math:: \Omega = 125.04452 - 1934.136261T + 0.0020708T^2 + T^3 / 450000

    The mean longitudes of the Sun, :math:`L`, and Moon, :math:`L'`, can be
    calculated in decimal degrees using the following:

    .. math:: L = 280.4665 + 36000.7698 * T

    .. math:: L' = 218.3165 + 481267.8813 * T

    The nutation in longitude, :math:`\Delta\psi`, and nutation in obliquity,
    :math:`\Delta\epsilon`, can then be calculated in decimal minutes. [Mee98c]_

    .. math:: \Delta\Psi = -17.20\sin\Omega - 1.32\sin 2L - 0.23\sin 2L' + 0.21\sin 2\Omega

    .. math:: \Delta\epsilon = 9.20\cos\Omega + 0.57\cos 2L + 0.10\cos 2L' - 0.09\cos 2\Omega

    References
    ----------
    .. [Mee98c] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 143–144. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> delta_psi, delta_epsilon = nutation_componenets(julian=julian)
    """

    T = (julian - 2451545) / 36525

    omega = 125.04452 - 1934.136261 * T + 0.0020708 * T ** 2 + T ** 3 / 450000
    omega = np.radians(omega)
    omega_2 = 2 * omega

    L = 2 * np.radians(280.4665 + 36000.7698 * T)
    LP = 2 * np.radians(218.3165 + 481267.8813 * T)

    P1 = -17.20 * np.sin(omega)
    P2 = -1.32 * np.sin(L)
    P3 = -0.23 * np.sin(LP)
    P4 = 0.21 * np.sin(omega_2)

    E1 = 9.20 * np.cos(omega)
    E2 = 0.57 * np.cos(L)
    E3 = 0.10 * np.cos(LP)
    E4 = -0.09 * np.cos(omega_2)

    delta_psi = P1 + P2 + P3 + P4
    delta_epsilon = E1 + E2 + E3 + E4

    return delta_psi, delta_epsilon

def mean_obliquity(julian: np.ndarray) -> np.ndarray:
    """Return the mean obliquity of the ecliptic in the J2000 epoch.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    np.ndarray
        Returns an array of shape (n,) containing the mean obliquity of the
        ecliptic in degrees and decimals.

    See Also
    --------
    apparent_obliquity :
        Return the apparent obliquity of the ecliptic in the J2000 epoch.

    Notes
    -----
    The time in Julian centeries since J2000, :math:`T`, can be calculated
    from the Julian day, :math:`JD`, from the following:

    .. math:: T = \\frac{JD - 2451545}{36525}

    If we define :math:`U = T / 100`, then the mean obliquity of the ecliptic.
    :math:`\epsilon_o` can be found by the following:

    .. math:: \epsilon_0 = 84381.448 - 4680.93 * U - 1.55 * U ** 2 +
        1999.25 * U ** 3 - 51.38 * U ** 4 - 249.67 * U ** 5 -
        39.05 * U ** 6 + 7.12 * U ** 7 + 27.87 * U ** 8 + 5.79 * U ** 9 +
        2.45 * U ** 10

    The algorithm used is only valid for 10000 years on either side of
    J2000. [Mee98d]_

    References
    ----------
    .. [Mee98d] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 147. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> epsilon_0 = mean_obliquity(julian=julian)
    """

    T = (julian - 2451545) / 36525
    U = T / 100

    epsilon_0 = 84381.448
    epsilon_0 -= 4680.93 * U
    epsilon_0 -= 1.55 * U ** 2
    epsilon_0 += 1999.25 * U ** 3
    epsilon_0 -= 51.38 * U ** 4
    epsilon_0 -= 249.67 * U ** 5
    epsilon_0 -= 39.05 * U ** 6
    epsilon_0 += 7.12 * U ** 7
    epsilon_0 += 27.87 * U ** 8
    epsilon_0 += 5.79 * U ** 9
    epsilon_0 += 2.45 * U ** 10

    epsilon_0 = epsilon_0 / 3600

    return epsilon_0

def apparent_obliquity(julian: np.ndarray) -> np.ndarray:
    """Return the apparent obliquity of the ecliptic in the J2000 epoch.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    np.ndarray
        Returns an array of shape (n,) containing the apparent obliquity of
        the ecliptic in degrees and decimals.

    See Also
    --------
    mean_obliquity :
        Return the mean obliquity of the ecliptic in the J2000 epoch.

    Notes
    -----
    The apparent obliquity of the ecliptic, :math:`\epsilon` can be calculated
    as :math:`\epsilon = \epsilon_0 + \Delta\epsilon` where :math:`\epsilon_0`
    is the mean obliquity of the ecliptic and :math:`\Delta\epsilon` is the
    nutation of obliquity. [Mee98e]_

    The due to the limitations of the algorithm used in the calculation of
    :math:`\epsilon_0` is only valid for 10000 years on either side of J2000.

    References
    ----------
    .. [Mee98e] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
       1998, pp. 147. isbn: 9780943396613.

    Examples
    --------
    >>> julian = np.array([2446895.5])
    >>> epsilon = apparent_obliquity(julian=julian)
    """

    epsilon_0 = mean_obliquity(julian)
    _, delta_epsilon = nutation_components(julian)
    epsilon = epsilon_0 + delta_epsilon / 3600

    return epsilon

def from_julian(julian: np.ndarray) -> Tuple:
    """Return year, month, day from Julian day data.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    Tuple
        Returns a tuple, `(year, month, day)`, where the items are NumPy
        arrays of shape (n,) containing the year, month, and decimal days of
        the input Julian data.

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
    ind_1 = np.where(month < 14)[0]
    ind_2 = np.where((month == 14) | (month == 15))[0]
    month[ind_1] = month[ind_1] - 1
    month[ind_2] = month[ind_2] - 13

    year = C
    ind_1 = np.where(month > 2)[0]
    ind_2 = np.where((month == 1) | (month == 2))[0]
    year[ind_1] = year[ind_1] - 4716
    year[ind_2] = year[ind_2] - 4715

    return year, month, day

def day_of_year(julian: np.ndarray) -> np.ndarray:
    """Return the day of the year.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing the day of the year.

    Notes
    -----
    The day of the year, :math:`N` can be calculated from the following:

    .. math:: N = INT\left(\frac{275M}{9}\right) - K\times\left(\frac{M+9}{12}\right) + D - 30

    where :math:`M` is the month number, :math:`D` is the day of the month,
    and :math:`K=1` if the year is a leap year otherwise :math:`K=2`. [Mee98g]_

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
    """Return the Equation of Time in degrees.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing the Equation of Time in decimal degrees.

    Notes
    -----
    The equation of time in radians can be calculated from the following:

    .. math:: E = y\sin 2L_0 - 2e\sin M + 4ey\sin M\cos 2L_0 - \frac{1}{2}y^2\sin 4L_0 - \frac{5}{4}e^2\sin 2M

    where :math:`y=\tan^2\frac{\epsilon}{2}`, :math:`\epsilon` is the apparent
    obliquity of the ecliptic, :math:`L_0` is the Sun's mean longitude,
    :math:`e` is the eccentricity of the Earth's orbit, and :math:`M` is the
    Sun's mean anomaly. [Mee98h]_

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

    T = (julian - 2451545) / 36525

    L0 = np.radians(280.46646 + 36000.76983 * T + 0.0003032 * T ** 2)
    M = np.radians(357.52911 + 35999.05029 * T - 0.0001537 * T ** 2)
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T ** 2

    epsilon = apparent_obliquity(julian=julian)
    y = np.tan(np.radians(epsilon / 2)) ** 2

    EOT = y * np.sin(2 * L0)
    EOT = EOT - 2 * e * np.sin(M)
    EOT = EOT + 4 * e * y * np.sin(M) * np.cos(2 * L0)
    EOT = EOT - 0.5 * y ** 2 * np.sin(4 * L0)
    EOT = EOT - 1.25 * e ** 2 * np.sin(2 * M)
    EOT = np.degrees(EOT)

    return EOT

def equation_of_equinoxes(julian: np.ndarray) -> np.ndarray:
    """Return the equation of the equinoxes in arcseconds.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing the Equation of Equinoxes in decimal
        arcseconds.

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

    delta_psi, _ = nutation_components(julian)
    epsilon = apparent_obliquity(julian)
    EoE = delta_psi * np.cos(np.radians(epsilon)) / 15

    return EoE


def sun_right_ascension(julian: np.ndarray) -> np.ndarray:
    """Return the right ascension of the mean sun position.

    This method calculates the Sun's right ascension to an accuracy of 0.01 of
    a degree.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing Julian times in decimal days.

    Returns
    -------
    np.ndarray
        Array of shape (n,) containing the right ascension of the mean sun in
        degrees and decimals.

    Notes
    -----
    To calculate the right ascension of the Sun, we must first calculate the
    following dependencies:

    .. math:: L0 = 280^\circ.46646 + 36000^\circ.76983t_U + 0^\circ.0003032t_U^2
    .. math:: M = 357^\circ.52911 + 35999^\circ.05029t_U - 0^\circ.0001537t_U^2

    .. math:: C = (1^\circ.914602 - 0^\circ.004817t_U - 0^\circ.000014t_U^2)\sin(m)
                + (0^\circ.019993 - 0^\circ.000101t_U)\sin(2m)
                + 0^\circ.000289\sin(3m)

    \odot = L0 + C

    where :math:`t_U = \frac{JD-2451545.0}{36525}` and `JD` is the Julian day.
    We can then calculate the right ascension, :math:`\alpha`, as

    .. math:: \tan\alpha = \frac{\cos\epsilon\sin\odot}{\cos\odot}

    where :math:`\epsilon` is the mean obliquity of the ecliptic. [Mee98j]_

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

    t_U = (julian - 2451545) / 36525
    t_U2 = t_U * t_U

    L0 = 280.46646 + 36000.76983 * t_U + 0.0003032 * t_U2
    M = 357.52911 + 35999.05029 * t_U - 0.0001537 * t_U2

    m = np.deg2rad(M)
    C = (1.914602 - 0.004817 * t_U - 0.000014 * t_U2) * np.sin(m)
    C = C + (0.019993 - 0.000101 * t_U) * np.sin(2 * m)
    C = C + 0.000289 * np.sin(3 * m)

    dot = np.radians(L0 + C)

    epsilon = mean_obliquity(julian=julian)
    epsilon = np.radians(epsilon)

    alpha = np.arctan2(np.cos(epsilon) * np.sin(dot), np.cos(dot))
    alpha = np.rad2deg(alpha) % 360

    return alpha
