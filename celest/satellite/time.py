

from celest.satellite._astronomical_quantities import (
    equation_of_time, from_julian, sun_right_ascension, equation_of_equinoxes
)
from datetime import datetime, timedelta
import numpy as np


class Time:
    """Time(julian, offset=0)

    Time transformations.

    `julian + offset` is the Julian time in J2000 epoch which can be converted
    into different time representations useful for astronomical calculations.

    Parameters
    ----------
    julian : array_like
        1-D array containing time data in Julian days.

        If time data is not in the J2000 epoch, an offset must be pased in to
        be added to the julian times.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch, default to zero.

    Examples
    --------
    Construct a `Time` object using J2000 input times:

    >>> time = [2460462.50, 2460462.91, 2460463.02]
    >>> t = Time(julian=time, offset=0)
    >>> gmst = t.gmst()

    For modified Julian input, we add an offset `2430000`:

    >>> time = [30462.50, 30462.91, 30463.02]
    >>> t = Time(julian=time, offset=2430000)

    We can then convert to various time representations:

    >>> gmst = t.gmst()
    >>> last = t.last(longitude=145)
    """

    def __init__(self, julian, offset=0) -> None:

        time = np.array(julian)

        if time.ndim == 0:
            time = time.reshape((1))
        elif time.ndim != 1:
            raise ValueError("Time array must be 1-D.")

        self._julian = time + offset
        self._length = self._julian.size

    def __len__(self):

        return self._length

    def true_solar_time(self, longitude) -> np.ndarray:
        """Return true solar time in decimal hours.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        np.ndarray
            1-D array containing true solar time in decimal hours.

        See Also
        --------
        mean_solar_time :
            Return mean solar time in decimal hours.

        Notes
        -----
        The true solar time, :math:`TT_s` can be calculated by converting
        Julian times, :math:`JD`, into UTC as follows:

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
        Calculate the true solar time for a constant longitude:

        >>> Time(julian=[2456293.54167]).true_solar_time(longitude=147.46)
        np.array([10.773109])
        """

        EoT = equation_of_time(julian=self._julian)

        idx = np.where(self._julian % 1 < 0.5)

        UTC = 24 * (self._julian % 1) - 12
        UTC[idx] = UTC[idx] + 12
        TST = (UTC + (np.array(longitude) + EoT) / 15) % 24

        return TST

    def mean_solar_time(self, longitude) -> np.ndarray:
        """Return mean solar time in decimal hours.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        np.ndarray
            1-D array containing mean solar time in decimal hours.

        See Also
        --------
        true_solar_time :
            Return true solar time in decimal hours.

        Notes
        -----
        The mean solar time, :math:`MT_s` can be calculated by converting
        Julian times, :math:`JD`, into UTC as follows:

        .. math:: UTC = 24\left(JD\%1\\right)^h + \\alpha^h

        The mean solar time can then be calculated using the Equation of Time,
        :math:`EoT`, and the local meridians longitude, :math:`\phi` as
        follows:

        .. math:: MT_s = (UTC^h + \phi^\circ/15)\%24

        Examples
        --------
        Calculate the mean solar time for a constant longitude:

        >>> Time(julian=2456293.54167).mean_solar_time(longitude=147.46)
        np.array([10.830667])
        """

        idx = np.where(self._julian % 1 < 0.5)

        UTC = 24 * (self._julian % 1) - 12
        UTC[idx] = UTC[idx] + 12
        MST = (UTC + np.array(longitude) / 15) % 24

        return MST

    def true_hour_angle(self, longitude) -> np.ndarray:
        """Return true hour angle in decimal hours.

        The true hour angle is the angle between the Sun's apparent position at
        a given time and its position at local solar noon. The value falls
        between `0` and `360`, measured in increasing degrees past local solar
        noon.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        np.ndarray
            1-D array containing true solar hour angles in decimal hours.

        See Also
        --------
        mean_hour_angle :
            Return mean hour angle in decimal hours.

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
        Calculate the true hour angle at different longitudes:

        >>> julian = [2455368.75, 2459450.85]
        >>> lon = [-105, -118.24]
        >>> Time(julian=julian).true_hour_angle(longitude=lon)
        np.array([10.97157662,
                  12.47807655])
        """

        TST = self.true_solar_time(longitude=longitude)
        HRA = (TST - 12) % 24

        return HRA

    def mean_hour_angle(self, longitude) -> np.ndarray:
        """Return mean hour angle in decimal hours.

        The mean hour angle is the angle between the Sun's mean position at
        a given time and its position at local solar noon. The value falls
        between `0` and `360`, measured in increasing degrees past local solar
        noon.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        np.ndarray
            1-D containing mean solar hour angles in decimal hours.

        See Also
        --------
        true_hour_angle :
            Return true hour angle in decimal hours.

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
        Calculate the meann hour angle for a single longitude:

        >>> Time(julian=2456293.5416666665).mean_hour_angle(longitude=147.46)
        np.array([22.83066666])
        """

        MST = self.mean_solar_time(longitude=longitude)
        HRA = (MST - 12) % 24

        return HRA

    def ut1(self) -> np.ndarray:
        """Return the universal time (same as GMT) in decimal hours.

        Due to approximations in the mean solar time calculations, the DUT1
        time correction is not accounted for which will introduce an error of
        at most 1.8 seconds.

        Returns
        -------
        np.ndarray
            1-D array containing universal time in decimal hours.

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
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
        >>> Time(julian=julian).ut1()
        np.array([6.00000,
                  8.40000,
                  1.00000])
        """

        ut1 = self.mean_solar_time(longitude=0)

        return ut1

    def julian(self) -> np.ndarray:
        """Convenience method to return Julian times.

        Returns
        -------
        np.ndarray
            1-D array containing the Julian times in decimal days.

        Examples
        --------
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
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
            1-D array containing `datetime.datetime` objects.

        Examples
        --------
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
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
        """Return Greenwich Mean Sidereal Time in decimal hours.

        Returns
        -------
        np.ndarray
            1-D array containing Greenwich Mean Sidereal Time in hours and
            decimals.

        Notes
        -----
        The Greenwich Mean Sidereal Time is calculated using the methods
        described in Astronomical Algorithms by Jean Meeus. [Mee98a]_

        References
        ----------
        .. [Mee98a] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 87 - 88. isbn: 9780943396613.

        Examples
        --------
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
        >>> Time(julian=julian).gmst()
        np.array([23.955316,
                  6.589391,
                  7.723214])
        """

        T = (self._julian - 2451545) / 36525
        T2 = T * T
        T3 = T2 * T

        gmst = 280.46061837
        gmst = gmst + 360.98564736629 * (self._julian - 2451545)
        gmst = gmst + 0.000387933 * T2
        gmst = gmst - T3 / 38710000
        gmst = gmst % 360 / 15

        return gmst

    def lmst(self, longitude) -> np.ndarray:
        """Return Local Mean Sidereal Time in decimal hours.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        np.ndarray
            1-D array containing Local Mean Sidereal Time in decimal hours.

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
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
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
        """Return Greenwich Apparent Sidereal Time in decimal hours.

        Returns
        -------
        np.ndarray
            1-D array containing Greenwich Apparent Sidereal Time in hours and
            degrees.

        Notes
        -----
        The difference between apparent and mean solar position is the equation
        of time. Thus, the Greenwich Apparent Sidereal Time can be calculated
        as follows:

        .. math:: GAST^h = GMST^h + EoT^h

        where :math:`EoT` is the equation of time in decimal hours.

        Examples
        --------
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
        >>> Time(julian=julian).gast()
        np.array([23.955596,
                  6.589146,
                  7.723465])
        """

        gmst = self.gmst()
        EoE = equation_of_equinoxes(self._julian) / 3600

        gast = gmst + EoE

        return gast

    def last(self, longitude) -> np.ndarray:
        """Return Local Apparent Sidereal Time in hours and degrees.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        np.ndarray
            1-D array containing Local Apparent Sidereal Time in hours and
            degrees.

        Notes
        -----
        The Local Apparent Sidereal Time is calculated using the methods
        described in Space-Time Reference Systems by M. Soffel and R.
        Langhans. [SL13f]_

        References
        ----------
        .. [SL13f] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 205.

        Examples
        --------
        >>> julian = [2455368.75, 2459450.85, 2456293.5416666665]
        >>> Time(julian=julian).last(longitude=150)
        np.array([9.98447567,
                  16.6286799,
                  17.78148932])
        """

        ut1 = self.ut1()
        alpha_sun = sun_right_ascension(self._julian) / 15
        EoE = equation_of_equinoxes(self._julian) / 3600

        last = (ut1 - 12 + alpha_sun + EoE + np.array(longitude) / 15) % 24

        return last
