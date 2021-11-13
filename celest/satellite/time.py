"""Time representations.

The time module contains the `Time` class to allow users to input an array of
times in Julian days and convert to various time representations.
"""


from celest.core.decorators import set_module
from celest.satellite._astronomical_quantities import (
    equation_of_time, from_julian, sun_right_ascension, equation_of_equinoxes
)
from datetime import datetime, timedelta
import numpy as np


@set_module('celest.satellite')
class Time(object):
    """Time representations.

    The `Time` class allows a user to convert an array of Julian day times in
    the J2000 epoch into various time representations.

    Parameters
    ----------
    julian : np.ndarray
        Array of shape (n,) containing time data in Julian days.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch.

    Methods
    -------
    true_solar_time(longitude)
        Return the true solar time (TTs) in hours and decimals.
    mean_solar_time(longitude)
        Return the mean solar time (MTs) in hours and decimals.
    true_hour_angle(longitude)
        Return the true hour angle in hours and decimals.
    mean_hour_angle(longitude)
        Return the mean hour angle in hours and decimals.
    ut1()
        Return the universal time (same as GMT) in hours and decimals.
    julian()
        Return Julian times.
    datetime()
        Return `datetime.datetime` object array.
    gmst(longitude)
        Return Greenwich Mean Sidereal Time in hours and decimals.
    lmst(longitude)
        Return Local Mean Sidereal Time in hours and decimals.
    gast(longitude)
        Return Greenwich Apparent Sidereal Time in hours and decimals.
    last(longitude)
        Return Local Apparent Sidereal Time in hours and degrees.
    """

    def __init__(self, julian: np.ndarray, offset: float=0) -> None:
        """Initialize attributes."""

        self._julian = julian + offset
        self._length = self._julian.size

    def __len__(self):
        """Return data size."""

        return self._length

    def true_solar_time(self, longitude: np.ndarray) -> np.ndarray:
        """Return the true solar time (TTs) in hours and decimals.

        Parameters
        ----------
        longitude : np.ndarray
            Array of shape (n,) containing longitude in decimal degrees.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing true solar time in hours and
            decimals.

        See Also
        --------
        mean_solar_time :
            Return the mean solar time (MTs) in hours and decimals.

        Notes
        -----
        The true solar time, :math:`TT_s` can be calculated by first converting
        Julian data, :math:`JD`, into UTC data by the following:

        .. math:: UTC = 24\left(JD\%1\\right)^h + \\alpha^h

        The true solar time can then be calculated using the Equation of Time,
        :math:`EoT`, definition and the mean solar time, :math:`MT_s`. [SL13a]_

        .. math:: TT_s = MT_s + EoT

        .. math:: TT_s = (UTC^h + (EoT^\circ + \phi^\circ)/15)\%24

        References
        ----------
        .. [SL13a] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julian = np.array([2456293.54167])
        >>> longitude = np.array([147.46])
        >>> Time(julian=julian).true_solar_time(longitude=longitude)
        np.array([10.773109])
        """

        julian = self._julian

        EoT = equation_of_time(julian=julian)

        alpha = np.full((julian.size,), -12)
        alpha[np.where(julian % 1 < 0.5)] = 12

        UTC = 24 * (julian % 1) + alpha
        TST = (UTC + (longitude + EoT) / 15) % 24

        return TST

    def mean_solar_time(self, longitude: np.ndarray) -> np.ndarray:
        """Return the mean solar time (MTs) in hours and decimals.

        Parameters
        ----------
        longitude : np.ndarray
            Array of shape (n,) containing longitude in degrees and decimals

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing mean solar time in hours and
            decimals.

        See Also
        --------
        true_solar_time :
            Return the true solar time (TTs) in hours and decimals.

        Notes
        -----
        The mean solar time, :math:`MT_s` can be calculated by first converting
        Julian data, :math:`JD`, into UTC data by the following:

        .. math:: UTC = 24\left(JD\%1\\right)^h + \\alpha^h

        The mean solar time can then be calculated using the Equation of Time,
        :math:`EoT`, and the local meridians longitude, :math:`\phi` as
        follows:

        .. math:: MT_s = (UTC^h + \phi^\circ/15)\%24

        Examples
        --------
        >>> julian = np.array([2456293.54167])
        >>> longitude = np.array([147.46])
        >>> Time(julian=julian).mean_solar_time(longitude=longitude)
        np.array([10.830667])
        """

        julian = self._julian

        alpha = np.full((julian.size,), -12)
        alpha[np.where(julian % 1 < 0.5)] = 12

        UTC = 24 * (julian % 1) + alpha
        MST = (UTC + longitude / 15) % 24

        return MST

    def true_hour_angle(self, longitude: np.ndarray) -> np.ndarray:
        """Return the true hour angle in hours and decimals.

        The true hour angle is the angle between the Sun's apparent position at
        the given times and its position at local solar noon. The value falls
        within the range :math:`[0, 360)` and is measured in increasing degrees
        past local solar noon.

        Parameters
        ----------
        longitude : np.ndarray
            Array of shape (n,) containing longitude in degrees and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing true solar hour angles in hours
            and decimals.

        See Also
        --------
        mean_hour_angle :
            Return the mean hour angle in hours and decimals.

        Notes
        -----
        The true hour angle is zero at local solar noon and increases with the
        true solar time by a rate of :math:`15^\circ` per hour. Thus, the true
        hour angle can be calculated from the following:

        .. math:: h_{Sun}^\circ = TT_s^\circ - 180^\circ

        where :math:`TT_s^\circ` is the true solar time in degrees. [SL13b]_

        References
        ----------
        .. [SL13b] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85])
        >>> lon = np.array([-105, -118.24])
        >>> Time(julian=julian).true_hour_angle(longitude=lon)
        np.array([10.97157662,
                  12.47807655])
        """

        TST = self.true_solar_time(longitude=longitude)
        HRA = (TST - 12) % 24

        return HRA

    def mean_hour_angle(self, longitude: np.ndarray) -> np.ndarray:
        """Return the mean hour angle in hours and decimals.

        The mean hour angle is the angle between the Sun's mean position at
        the given times and its position at local solar noon. The value falls
        within the range :math:`[0, 360)` and is measured in increasing degrees
        past local solar noon.

        Parameters
        ----------
        longitude : np.ndarray
            Array of shape (n,) containing longitude in degrees and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing mean solar hour angles in hours and
            decimals.

        See Also
        --------
        true_hour_angle :
            Return the true hour angle in hours and decimals.

        Notes
        -----
        The mean hour angle is zero at local solar noon and increases with the
        mean solar time by a rate of :math:`15^\circ` per hour. Thus, the mean
        hour angle can be calculated from the following:

        .. math:: h_{hSun}^\circ = MT_s^\circ - 180^\circ

        where :math:`MT_s^\circ` is the mean solar time in degrees. [SL13c]_

        References
        ----------
        .. [SL13c] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julian = np.array([2456293.5416666665])
        >>> longitude = np.array([147.46])
        >>> Time(julian=julian).mean_hour_angle(longitude=longitude)
        np.array([22.83066666])
        """

        MST = self.mean_solar_time(longitude=longitude)
        HRA = (MST - 12) % 24

        return HRA

    def ut1(self) -> np.ndarray:
        """Return the universal time (same as GMT) in hours and decimals.

        This method returns a universal time by calculating the mean solar time
        at the Greenwich prime meridian. Due to approximations in the mean
        solar time calculations, the DUT1 time correction is not accounted for
        which will introduce an error of at most 1.8 seconds.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing universal time in hours and
            decimals.

        Notes
        -----
        The universal time is equal to the mean solar time at the Greenwich
        meridian. [SL13d]_ It can be calculated using
        `Time.mean_solar_time(longitude=0)`.

        References
        ----------
        .. [SL13d] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).ut1()
        np.array([6.00000,
                  8.40000,
                  1.00000])
        """

        ut1 = self.mean_solar_time(longitude=0)

        return ut1

    def julian(self) -> np.ndarray:
        """Return Julian times.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing the Julian times in days and
            decimals.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).julian()
        np.array([2455368.75,
                  2459450.85,
                  2456293.5416666665])
        """

        return self._julian

    def datetime(self) -> np.ndarray:
        """Return `datetime.datetime` object array.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing `datetime.datetime` objects.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).datetime()
        np.array([datetime.datetime(2010, 6, 21, 6, 0)
                  datetime.datetime(2021, 8, 24, 8, 24, 0, 8)
                  datetime.datetime(2013, 1, 1, 0, 59, 59, 999987)])
        """

        julian = self._julian

        year, month, day = from_julian(julian)
        remainder = day % 1
        day = day.astype(int)

        datetime_arr = np.zeros(julian.shape).astype("O")
        for i in range(julian.size):
            y, m, d, r = year[i], month[i], day[i], remainder[i]
            dt = datetime(year=y, month=m, day=d) + timedelta(days=r)
            datetime_arr[i] = dt

        return datetime_arr

    def gmst(self) -> np.ndarray:
        """Return Greenwich Mean Sidereal Time in hours and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing Greenwich Mean Sidereal Time in
            hours and decimals.

        Notes
        -----
        The Greenwich Mean Sidereal Time can be calculated from the following
        equation:

        .. math:: GMST = \left(280.46 + 360.99(j - 2451545) + 0.000387933T^2 - \\frac{T^3}{38710000}\\right) % 360 / 15

        where

        .. math:: T = \\frac{JD - 2451545}{36525}

        and :math:`JD` is the current Julian day. [Mee98a]_

        References
        ----------
        .. [Mee98a] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 87 - 88. isbn: 9780943396613.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).gmst()
        np.array([23.955316,
                  6.589391,
                  7.723214])
        """

        julData = self._julian

        T = (julData - 2451545) / 36525
        T2 = T * T
        T3 = T2 * T

        gmst = 280.46061837
        gmst = gmst + 360.98564736629 * (julData - 2451545)
        gmst = gmst + 0.000387933 * T2
        gmst = gmst - T3 / 38710000
        gmst = gmst % 360 / 15

        return gmst

    def lmst(self, longitude: np.ndarray) -> np.ndarray:
        """Return Local Mean Sidereal Time in hours and decimals.

        Parameters
        ----------
        longitude : np.ndarray
            Array of shape (n,) containing longitude in degrees and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing Local Mean Sidereal Time in hours and
            decimals.

        Notes
        -----
        The Local Mean Sidereal Time can be calculated using the following
        formulation:

        .. math:: LMST^h = h^h_{mSun} + \\alpha^h_{mSun}

        where :math:`h^h_{mSun}` is the mean solar hour angle at the observer's
        longitude and :math:`\\alpha^h_{mSun}` is the right ascension of the
        mean Sun position. [SL13e]_

        References
        ----------
        .. [SL13e] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 205.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).lmst(longitude=150)
        np.array([9.98419514,
                  16.62892539,
                  17.78123885])
        """

        hour_angle = self.mean_hour_angle(longitude)
        alpha_sun = sun_right_ascension(self._julian) / 15

        lmst = (hour_angle + alpha_sun) % 24

        return lmst

    def gast(self) -> np.ndarray:
        """Return Greenwich Apparent Sidereal Time in hours and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing Greenwich Apparent Sidereal Time in
            hours and degrees.

        Notes
        -----
        The difference between apparent and mean solar position is the equation
        of time. Thus, the Greenwich Apparent Sidereal Time can be calculated
        as follows:

        .. math:: GAST^h = GMST^h + EoT^h

        where :math:`EoT` is the equation of time in hours and decimals.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).gast()
        np.array([23.955596,
                  6.589146,
                  7.723465])
        """

        gmst = self.gmst()
        EoE = equation_of_equinoxes(self._julian) / 3600

        gast = gmst + EoE

        return gast

    def last(self, longitude: np.ndarray) -> np.ndarray:
        """Return Local Apparent Sidereal Time in hours and degrees.

        Parameters
        ----------
        longitude : np.ndarray
            Array of shape (n,) containing the actual astronomical longitude
            in degrees and decimals.

        Returns
        -------
        np.ndarray
            Array of shape (n,) containing Local Apparent Sidereal Time in
            hours and degrees.

        Notes
        -----
        The Local Apparent Sidereal Time can be calculated using the following
        formulation:

        .. math:: LAST^h = UT1^h - 12^h + \\alpha^h_{mSun} + EoE^h + \Lambda^\circ/15

        where :math:`\\alpha^h_{mSun}` is the right ascension of the mean Sun
        position in hours, :math:`EoE^h` is the equation of equinoxes in hours,
        and :math:`\Lambda` is the observer's longitude in  degrees. [SL13f]_

        References
        ----------
        .. [SL13f] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 205.

        Examples
        --------
        >>> julian = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julian).last(longitude=150)
        np.array([9.98447567,
                  16.6286799,
                  17.78148932])
        """

        ut1 = self.ut1()
        alpha_sun = sun_right_ascension(self._julian) / 15
        EoE = equation_of_equinoxes(self._julian) / 3600

        last = (ut1 - 12 + alpha_sun + EoE + longitude / 15) % 24

        return last
