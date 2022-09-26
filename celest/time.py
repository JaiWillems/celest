

from celest.constants import (
    DAYS_IN_JULIAN_CENTURY,
    ERA_AT_JD2000_DEG
)
from celest.coordinates.astronomical_quantities import (
    equation_of_time,
    from_julian,
    _calculate_raw_elapsed_jd_century_powers,
    sun_right_ascension,
    equation_of_equinoxes
)
from celest.units.quantity import Quantity
from celest import units as u
from datetime import datetime, timedelta
import numpy as np


class Time:
    """Time(julian, offset=0)

    Time transformations.

    `julian + offset` is the Julian time in JD2000 epoch which can be converted
    into different time representations useful for astronomical calculations.

    Parameters
    ----------
    julian : array_like
        1-D array containing time data in Julian days.

        If time data is not in the JD2000 epoch, an offset must be passed in to
        be added to the julian times.
    offset : float, optional
        Offset to convert input time data to the JD2000 epoch, default to zero.

    Methods
    -------
    true_solar_time(longitude)
        Return true solar time.
    mean_solar_time(longitude)
        Return mean solar time.
    true_hour_angle(longitude)
        Return true hour angle.
    mean_hour_angle(longitude)
        Return mean hour angle.
    ut1()
        Return the universal time (same as GMT).
    datetime()
        Return `datetime.datetime` object array.
    gmst()
        Return Greenwhich Mean Sidereal Time.
    lmst(longitude)
        Return Local Mean Sidereal Time.
    gast()
        Return Greenwich Apparent Sidereal Time.
    last(longitude)
        Return Local Apparent Sidereal Time.

    Examples
    --------
    Construct a `Time` object using JD2000 input times:

    >>> julian = [2460462.50, 2460462.91, 2460463.02]
    >>> time = Time(julian=time, offset=0)
    >>> gmst = time.gmst()

    For modified Julian input, we add an offset `2430000`:

    >>> time = [30462.50, 30462.91, 30463.02]
    >>> time = Time(julian=time, offset=2430000)

    We can then convert to various time representations:

    >>> gmst = time.gmst()
    >>> last = time.last(longitude=145)
    """

    def __init__(self, julian: np.ndarray, offset: float=0) -> None:
        """Time transformations.

        `julian + offset` is the Julian time in JD2000 epoch which can be
        converted into different time representations useful for astronomical
        calculations.

        Parameters
        ----------
        julian : array_like
            1-D array containing time data in Julian days.

            If time data is not in the JD2000 epoch, an offset must be passed in
            to be added to the julian times.
        offset : float, optional
            Offset to convert input time data to the JD2000 epoch, default to
            zero.
        """

        time = np.array(julian)

        if time.ndim == 0:
            time = time.reshape((1,))
        elif time.ndim != 1:
            raise ValueError("Time array must be 1-D.")

        self._julian = time + offset
        self._length = self._julian.size

    def __len__(self):
        return self._length

    @property
    def julian(self) -> Quantity:
        """
        Return Julian times.
        """

        return Quantity(self._julian, u.jd2000)

    def true_solar_time(self, longitude: np.ndarray) -> Quantity:
        """Return true solar time.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        Quantity
            Quantity object containing true solar time.

        See Also
        --------
        mean_solar_time :
            Return mean solar time.

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
        Quantity(np.array([10.773109]), Unit("hourangle"))
        """

        eqn_of_time = equation_of_time(self.julian).to(u.deg).data
        utc_time = 24 * (self._julian % 1) - 12
        true_solar_time = (utc_time + (np.array(longitude) + eqn_of_time) / 15) % 24

        return Quantity(true_solar_time, u.hourangle)

    def mean_solar_time(self, longitude: np.ndarray) -> Quantity:
        """Return mean solar time.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        Quantity
            Quantity object containing mean solar time.

        See Also
        --------
        true_solar_time :
            Return true solar time.

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
        Quantity(np.array([10.830667]), Unit("hourangle"))
        """

        utc_time = 24 * (self._julian % 1) - 12
        mean_solar_time = (utc_time + np.array(longitude) / 15) % 24

        return Quantity(mean_solar_time, u.hourangle)

    def true_hour_angle(self, longitude) -> Quantity:
        """Return true hour angle.

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
        Quantity
            Quantity object containing true solar hour angles.

        See Also
        --------
        mean_hour_angle :
            Return mean hour angle.

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
        >>> longitude = [-105, -118.24]
        >>> Time(julian=julian).true_hour_angle(longitude=longitude)
        Quantity(np.array([10.97157662, 12.47807655]), Unit("hourangle"))
        """

        tha = (self.true_solar_time(longitude).to(u.hourangle) - 12) % 24
        return Quantity(tha, u.hourangle)

    def mean_hour_angle(self, longitude: np.ndarray) -> Quantity:
        """Return mean hour angle.

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
        Quantity
            Quantity object containing mean solar hour angles.

        See Also
        --------
        true_hour_angle :
            Return true hour angle.

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
        Quantity(np.array([22.83066666]), Unit("hourangle"))
        """

        mha = (self.mean_solar_time(longitude).to(u.hourangle) - 12) % 24
        return Quantity(mha, u.hourangle)

    def ut1(self) -> Quantity:
        """Return the universal time (same as GMT).

        Due to approximations in the mean solar time calculations, the DUT1
        time correction is not accounted for which will introduce an error of
        at most 1.8 seconds.

        Returns
        -------
        Quantity
            Quantity object containing universal time.

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
        Quantity(np.array([6.00000, 8.40000, 1.00000]), Unit("hourangle"))
        """

        return self.mean_solar_time(0)

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

        year, month, day = from_julian(Quantity(self._julian, u.jd2000))
        full_days, fractional_day = day.astype(int), day % 1

        datetime_list = []
        for y, m, d, f in zip(year, month, full_days, fractional_day):
            integer_date = datetime(year=y, month=m, day=d)
            extra_date = timedelta(days=f)
            datetime_list.append(integer_date + extra_date)

        return np.array(datetime_list)

    def gmst(self) -> Quantity:
        """Return Greenwich Mean Sidereal Time.

        Returns
        -------
        Quantity
            Quantity object containing Greenwich Mean Sidereal Time.

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
        Quantity(np.array([23.955316, 6.589391, 7.723214]), Unit("hourangle"))
        """

        T1, T2, T3 = _calculate_raw_elapsed_jd_century_powers(self.julian, 3)

        gmst = ERA_AT_JD2000_DEG + 360.98564736629 * T1 *\
            DAYS_IN_JULIAN_CENTURY + 0.000387933 * T2 - T3 / 38710000

        return Quantity(gmst % 360 / 15, u.hourangle)

    def lmst(self, longitude: np.ndarray) -> Quantity:
        """Return Local Mean Sidereal Time.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        Quantity
            Quantity object containing Local Mean Sidereal Time.

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
        Quantity(np.array([9.98419514, 16.62892539, 17.78123885]), Unit("hourangle"))
        """

        hour_angle = self.mean_hour_angle(longitude).to(u.hourangle)
        right_ascension = sun_right_ascension(self.julian).to(u.hourangle)

        return Quantity((hour_angle + right_ascension) % 24, u.hourangle)

    def gast(self) -> Quantity:
        """Return Greenwich Apparent Sidereal Time.

        Returns
        -------
        Quantity
            Quantity object containing Greenwich Apparent Sidereal Time.

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

        gmst = self.gmst().to(u.hourangle)
        eqn_of_equinoxes = equation_of_equinoxes(self.julian).to(u.deg)

        return Quantity(gmst + eqn_of_equinoxes, u.hourangle)

    def last(self, longitude: np.ndarray) -> Quantity:
        """Return Local Apparent Sidereal Time.

        Parameters
        ----------
        longitude : array_like
            1-D array containing longitude in decimal degrees. If a scalar
            value is given, it will be applied for all times.

        Returns
        -------
        Quantity
            Quantity object containing Local Apparent Sidereal Time.

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

        ut1 = self.ut1().to(u.hourangle)
        right_ascension_of_sun = sun_right_ascension(self.julian).to(u.hourangle)
        eqn_of_equinoxes = equation_of_equinoxes(self.julian).to(u.deg)

        last = ut1 - 12 + right_ascension_of_sun + eqn_of_equinoxes +\
            np.array(longitude) / 15

        return Quantity(last % 24, u.hourangle)
