"""Time representations.

The time module contains the `Time` class to allow users to input an array of
times in julian days and convert to various time representations.
"""


from celest.core.decorators import set_module
from celest.satellite import AstronomicalQuantities, Interpolation
from datetime import datetime, timedelta
from typing import Any, Dict
import numpy as np


@set_module('celest.satellite')
class Time(AstronomicalQuantities, Interpolation):
    """Time representations.
    
    The `Time` class allows a user to convert an array of times in julian days
    into various time representations.

    Parameters
    ----------
    julian : np.array
        Array of shape (n,) containing time data in julian days.
    offset : float, optional
        Offset to convert input time data to the J2000 epoch.
    factor : int, optional
        Interpolate the inputted time data by `factor` times.
    
    Methods
    -------
    true_solar_time(longitude, **kwargs)
        Return the true solar time (TTs) in hours and decimals.
    mean_solar_time(longitude, **kwargs)
        Return the mean solar time (MTs) in hours and decimals.
    true_hour_angle(longitude, **kwargs)
        Return the true hour angle in hours and decimals.
    mean_hour_angle(longitude, **kwargs)
        Return the mean hour angle in hours and decimals.
    sun_right_ascension()
        Return the right ascension of the mean sun position.
    UT1(**kwargs)
        Return the universal time (same as GMT) in hours and decimals.
    julian(**kwargs)
        Return julian times.
    datetime(**kwargs)
        Return `datetime.datetime` object array.
    GMST(longitude, **kwargs)
        Return Greenwhich Mean Sidereal Time in hours and decimals.
    LMST(longitude, **kwargs)
        Return Local Mean Sidereal Time in hours and decimals.
    GAST(longitude, **kwargs)
        Return Greenwhich Apparent Sidereal Time in hours and decimals.
    LAST(longitude, **kwargs)
        Return Local Apparent Sidereal Time in hours and degrees.
    """

    def __init__(self, julian: np.array, offset: float=0, factor: int=0) -> None:
        """Initialize attributes."""

        if factor > 0:
            julian = self._interp(data=julian, factor=factor)

        self._julian = julian + offset
        self._length = self._julian.size
    
    def __len__(self):
        """Return data size."""

        return self._length

    def true_solar_time(self, longitude: np.array, **kwargs:
                        Dict[str, Any]) -> np.array:
        """Return the true solar time (TTs) in hours and decimals.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
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

        .. math:: UTC = 24\left(JD\%1\right)^h + \alpha^h

        The true solar time can then be calculated using the Equation of Time,
        :math:`EoT`, definition and the mean solar time, :math:`MT_s`.

        .. math:: TT_s = MT_s + EoT [1]_

        .. math:: TT_s = (UTC^h + (EoT^\circ + \phi^\circ)/15)\%24

        References
        ----------
        .. [1] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julData = np.array([2456293.54167])
        >>> longitude = np.array([147.46])
        >>> Time(julian=julData).true_solar_time(longitude=longitude)
        np.array([10.773109])
        """
        
        EoT = self.equation_of_time(julData=self._julian)

        julian = self._julian

        alpha = np.full((julian.size,), -12)
        alpha[np.where(julian % 1 < 0.5)] = 12

        UTC = 24 * (julian % 1) + alpha
        TST = (UTC + (longitude + EoT) / 15) % 24

        if kwargs:
            TST = self._interp(TST, **kwargs)

        return TST

    def mean_solar_time(self, longitude: np.array, **kwargs:
                        Dict[str, Any]) -> np.array:
        """Return the mean solar time (MTs) in hours and decimals.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in degrees and decimals
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
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

        .. math:: UTC = 24\left(JD\%1\right)^h + \alpha^h

        The mean solar time can then be calculated using the Equation of Time,
        :math:`EoT`, and the local meridians longitude, :math:`\phi` as
        follows:

        .. math:: MT_s = (UTC^h + \phi^\circ/15)\%24

        Examples
        --------
        >>> julData = np.array([2456293.54167])
        >>> longitude = np.array([147.46])
        >>> Time(julian=julData).mean_solar_time(longitude=longitude)
        np.array([10.830667])
        """

        julian = self._julian

        alpha = np.full((julian.size,), -12)
        alpha[np.where(julian % 1 < 0.5)] = 12

        UTC = 24 * (julian % 1) + alpha
        MST = (UTC + longitude / 15) % 24

        if kwargs:
            MST = self._interp(MST, **kwargs)

        return MST
    
    def true_hour_angle(self, longitude: np.array, **kwargs:
                        Dict[str, Any]) -> np.array:
        """Return the true hour angle in hours and decimals.

        The true hour angle is the angle between the Sun's apparent position at
        the given times and its position at local solar noon. The value falls
        within the range :math:`[0, 360)` and is measured in increasing degrees
        past local solar noon.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in degrees and decimals.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
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

        where :math:`TT_s^\circ` is the true solar time in degrees.[1]_

        References
        ----------
        .. [1] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85])
        >>> lon = np.array([-105, -118.24])
        >>> Time(julian=julData).true_hour_angle(longitude=lon)
        np.array([10.97157662,
                  12.47807655])
        """

        TST = self.true_solar_time(longitude=longitude)
        HRA = (TST - 12) % 24

        if kwargs:
            HRA = self._interp(HRA, **kwargs)

        return HRA

    def mean_hour_angle(self, longitude: np.array, **kwargs:
                        Dict[str, Any]) -> np.array:
        """Return the mean hour angle in hours and decimals.

        The mean hour angle is the angle between the Sun's mean position at
        the given times and its position at local solar noon. The value falls
        within the range :math:`[0, 360)` and is measured in increasing degrees
        past local solar noon.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in degrees and decimals.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
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

        where :math:`MT_s^\circ` is the mean solar time in degrees.[1]_

        References
        ----------
        .. [1] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julData = np.array([2456293.5416666665])
        >>> longitude = np.array([147.46])
        >>> Time(julian=julData).mean_hour_angle(longitude=longitude)
        np.array([22.83066666])
        """

        MST = self.mean_solar_time(longitude=longitude)
        HRA = (MST - 12) % 24

        if kwargs:
            HRA = self._interp(HRA, **kwargs)

        return HRA
    
    def sun_right_ascension(self) -> np.array:
        """Return the right ascension of the mean sun position.

        This method calculates the Sun's right ascension to an accuracy of 0.01
        of a degree.

        Returns
        -------
        np.array
            Array of shape (n,) containing the right ascension of the mean sun
            in degrees and decimals.
        
        Notes
        -----
        To calculate the right ascension of the Sun, we must first calculate
        the following dependencies:

        .. math:: L0 = 280^\circ.46646 + 36000^\circ.76983t_U + 0^\circ.0003032t_U^2
        .. math:: M = 357^\circ.52911 + 35999^\circ.05029t_U - 0^\circ.0001537t_U^2

        .. math:: C = (1^\circ.914602 - 0^\circ.004817t_U - 0^\circ.000014t_U^2)\sin(m)
                    + (0^\circ.019993 - 0^\circ.000101t_U)\sin(2m)
                    + 0^\circ.000289\sin(3m)

        \odot = L0 + C

        where :math:`t_U = \frac{JD-2451545.0}{36525}` and `JD` is the Julian
        day. We can then calculate the right ascension, :math:`\alpha`, as

        .. math:: \tan\alpha = \frac{\cos\epsilon\sin\odot}{\cos\odot}

        where :math:`\epsilon` is the mean obliquity of the ecliptic.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 163 - 165. isbn: 9780943396613.

        Examples
        --------
        >>> julData = np.array([2448908.5])
        >>> Time(julian=julData).sun_right_ascension()
        np.array([198.38166])
        """

        t_U = (self._julian - 2451545) / 36525
        t_U2 = t_U * t_U

        L0 = 280.46646 + 36000.76983 * t_U + 0.0003032 * t_U2
        M = 357.52911 + 35999.05029 * t_U - 0.0001537 * t_U2

        m = np.deg2rad(M)
        C = (1.914602 - 0.004817 * t_U - 0.000014 * t_U2) * np.sin(m)
        C = C + (0.019993 - 0.000101 * t_U) * np.sin(2 * m)
        C = C + 0.000289 * np.sin(3 * m)

        dot = np.radians(L0 + C)

        epsilon = self.mean_obliquity(julData=self._julian)
        epsilon = np.radians(epsilon)

        alpha = np.arctan2(np.cos(epsilon) * np.sin(dot), np.cos(dot))
        alpha = np.rad2deg(alpha) % 360

        return alpha

    def UT1(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the universal time (same as GMT) in hours and decimals.

        This method returns a universal time by calculating the mean solar time
        at the Greenwich prime meridian. Due to approcimations in the mean
        solar time calculations the DUT1 time correction is not accounted for
        which will introduce an error of at most 1.8 seconds.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing universal time in hours and
            decimals.

        Notes
        -----
        The universal time is equal to the mean solar time at the Greenwich
        meridian.[1]_ It can be calculated using
        `Time.mean_solar_time(longitude=0)`.

        References
        ----------
        .. [1] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 203.

        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).UT1()
        np.array([6.00000,
                  8.40000,
                  1.00000])
        """

        UT1 = self.mean_solar_time(longitude=0)

        if kwargs:
            UT1 = self._interp(UT1, **kwargs)

        return UT1
    
    def julian(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return julian times.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing the julian times in days and
            decimals.
        
        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).julian()
        np.array([2455368.75,
                  2459450.85,
                  2456293.5416666665])
        """

        julian = self._julian

        if kwargs:
            julian = self._interp(julian, **kwargs)

        return julian

    def datetime(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return `datetime.datetime` object array.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing `datetime.datetime` objects.
        
        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).datetime()
        np.array([datetime.datetime(2010, 6, 21, 6, 0)
                  datetime.datetime(2021, 8, 24, 8, 24, 0, 8)
                  datetime.datetime(2013, 1, 1, 0, 59, 59, 999987)])
        """

        jul_data = self._julian

        if kwargs:
            jul_data = self._interp(jul_data, **kwargs)

        year, month, day = self.from_julian(jul_data)
        remainder = day % 1
        day = day.astype(int)

        datetime_arr = np.zeros(jul_data.shape).astype("O")
        for i in range(jul_data.size):
            y, m, d, r = year[i], month[i], day[i], remainder[i]
            dt = datetime(year=y, month=m, day=d) + timedelta(days=r)
            datetime_arr[i] = dt

        return datetime_arr
    
    def GMST(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return Greenwhich Mean Sidereal Time in hours and decimals.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing Greenwhich Mean Sideral Time in
            hours and decimals.
        
        Notes
        -----
        The Greenwhich Mean Sidereal Time can be calculated from the following
        equation:

        .. math:: GMST = \left(280.46061837 + 360.98564736629(j - 2451545) +
            0.000387933T^2 - \\frac{T^3}{38710000}\\right) % 360 / 15

        where

        .. math:: T = \\frac{JD - 2451545}{36525}

        and :math:`JD` is the current Julian day.[1]_

        References
        ----------
        .. [1] Jean Meeus. Astronomical algorithms. 2nd ed. Willmann-Bell,
           1998, pp. 87 - 88. isbn: 9780943396613.

        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).GMST()
        np.array([23.955316,
                  6.589391,
                  7.723214])
        """

        julData = self._julian

        T = (julData - 2451545) / 36525
        T2 = T * T
        T3 = T2 * T

        GMST = 280.46061837
        GMST = GMST + 360.98564736629 * (julData - 2451545)
        GMST = GMST + 0.000387933 * T2
        GMST = GMST - T3 / 38710000
        GMST = GMST % 360 / 15

        if kwargs:
            GMST = self._interp(GMST, **kwargs)
        
        return GMST
    
    def LMST(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return Local Mean Sidereal Time in hours and decimals.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in degrees and decimals.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing Local Mean Sideral Time in hours and
            decimals.
        
        Notes
        -----
        The Local Mean Sidereal Time can be calculated using the following
        formulation:

        .. math:: LMST^h = h^h_{mSun} + \alpha^h_{mSun}

        where :math:`h^h_{mSun}` is the mean solar hour angle at the observer's
        longitude and :math:`\alpha^h_{mSun}` is the right ascension of the
        mean Sun position.[1]_

        References
        ----------
        .. [1] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 205.
        
        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).LMST(longitude=150)
        np.array([9.98419514,
                  16.62892539,
                  17.78123885])
        """

        hour_angle = self.mean_hour_angle(longitude)
        alpha_sun = self.sun_right_ascension() / 15

        LMST = (hour_angle + alpha_sun) % 24

        if kwargs:
            LMST = self._interp(LMST, **kwargs)
        
        return LMST
    
    def GAST(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return Greenwhich Apparent Sidereal Time in hours and decimals.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing Greenwhich Apparent Sideral Time in
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
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).GAST()
        np.array([23.955596,
                  6.589146,
                  7.723465])
        """

        GMST = self.GMST()
        EoE = self.equation_of_equinoxes(self._julian) / 3600

        GAST = GMST + EoE

        if kwargs:
            GAST = self._interp(GAST, **kwargs)

        return GAST

    def LAST(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return Local Apparent Sidereal Time in hours and degrees.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing the actual astronomical longitude
            in degrees and decimals.
        kwargs : Dict, optional
            Optional keyword arguments passed into the
            `Interpolation()._interp()` method.

        Returns
        -------
        np.array
            Array of shape (n,) containing Local Apparent Sideral Time in
            hours and degrees.
        
        Notes
        -----
        The Local Apparent Sidereal Time can be calculated using the following
        formulation:

        .. math:: LAST^h = UT1^h - 12^h + \alpha^h_{mSun} + EoE^h + \Lambda^\circ/15

        where :math:`\alpha^h_{mSun}` is the right ascension of the mean Sun
        position in hours, :math:`EoE^h` is the equation of equinoxes in hours,
        and :math:`\Lambda` is the observers longitude in  degrees.[1]_
        
        References
        ----------
        .. [1] M. Soffel and R. Langhans. Space-Time Reference Systems.
           Astronomy and Astrophysics Library. Springer-Verlag, 2013, p. 205.
        
        Examples
        --------
        >>> julData = np.array([2455368.75, 2459450.85, 2456293.5416666665])
        >>> Time(julian=julData).LAST(longitude=150)
        np.array([9.98447567,
                  16.6286799,
                  17.78148932])
        """

        UT1 = self.UT1()
        alpha_sun = self.sun_right_ascension() / 15
        EoE = self.equation_of_equinoxes(self._julian) / 3600

        LAST = (UT1 - 12 + alpha_sun + EoE + longitude / 15) % 24

        if kwargs:
            LAST = self._interp(LAST, **kwargs)

        return LAST
