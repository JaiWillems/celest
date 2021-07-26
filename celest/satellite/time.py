"""Time representations.

The time module contains the `Time` class to allow users to input an array of
times in julian days and convert to various time representations.
"""


from celest.astronomy import Sun
from celest.core.decorators import set_module
from celest.satellite import AstronomicalQuantities, Interpolation
from typing import Any, Dict
import julian
import numpy as np


@set_module('celest.satellite')
class Time(AstronomicalQuantities):
    """Time representations.
    
    The `Time` class allows a user to input an array of times in julian days
    and convert to various time representations.

    Parameters
    ----------
    julian : np.array
        Array of shape (n,) containing julian dates.
    offset : float, optional
        Offset to convert inputted julian dates to the J2000 epoch.
    factor : int, optional
        Interpolate the inputted time data by `factor` times.
    
    Attributes
    ----------
    julian : np.array
        Array of julian dates interpolated by `factor` times.
    length : int
        Length of the `julian` attribute.
    
    Methods
    -------
    true_hour_angle(longitude, **kwargs)
        Return the hour angle using the true sun position.
    mean_hour_angle(longitude, **kwargs)
        Return the hour angle of the mean sun position.
    true_solar_time(longitude, **kwargs)
        Return the true solar time (TTs).
    mean_solar_time(longitude, **kwargs)
        Return the mean solar time (MTs, same as LMT).
    sun_right_ascension()
        Return the right ascension of the mean sun position.
    UT1(**kwargs)
        Return the universal time (same as GMT).
    julian(**kwargs)
        Return julian times.
    datetime_UTC(**kwargs)
        Return `datetime.datetime` UTC strings.
    GMST(longitude, **kwargs)
        Return Greenwhich Mean Sidereal Time.
    LMST(longitude, **kwargs)
        Return Local Mean Sidereal Time.
    """

    def __init__(self, julian: np.array, offset: float=0, factor: int=0) -> None:
        """Initialize attributes."""

        self._equation_of_time = None
        self._sun_right_ascension = None
        self._UT1 = None
        self._equation_of_equinoxes = None
        self._GMST = None
        self._GAST = None

        self._interp = Interpolation()

        if factor > 0:
            julian = self._interp(data=julian, factor=factor)

        self._julian = julian + offset
        self.length = self._julian.size

    def true_hour_angle(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return the hour angle using the true sun position.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing true solar hour angles in decimal
            degrees.
        
        See Also
        --------
        mean_hour_angle
        """

        theta = np.radians(longitude)

        sun_position = Sun().position(self).ITRS(self)
        obs_position = np.full((self.length, 3), [np.cos(theta), np.sin(theta), 0])
        arb1 = np.full((self.length, 3), [0, 0, 100])
        arb2 = np.full((self.length, 3), [0, 0, -100])

        norm_a = np.cross(sun_position - arb2, sun_position - arb1)
        norm_b = np.cross(obs_position - arb2, obs_position - arb1)

        na_dot_nb = np.diag(np.matmul(norm_a, norm_b.T))
        denom = np.linalg.norm(norm_a, axis=1) * np.linalg.norm(norm_b, axis=1)
        HRA = np.degrees(np.arccos(na_dot_nb / denom))

        sun_dot_nb = np.sum(sun_position * norm_b, axis=1)
        denom = np.linalg.norm(norm_b, axis=1) ** 2
        scale_arr = np.repeat((sun_dot_nb / denom).reshape((-1, 1)), 3, axis=1)
        proj_v4_on_nb = scale_arr * norm_b

        test = np.sum(proj_v4_on_nb * norm_b, axis=1)

        pos_ind = np.where(test >= 0)[0]

        HRA[pos_ind] = -HRA[pos_ind]

        if kwargs:
            HRA = self._interp(HRA, **kwargs)

        return HRA

    def mean_hour_angle(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return the hour angle of the mean sun position.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing mean solar hour angles in decimal
            degrees.
        
        See Also
        --------
        true_hour_angle

        Notes
        -----
        The mean solar hour angle, :math:`h_{mSun}`, is found from the
        following:

        .. math:: h_{mSun} = h_{Sun} - Eq.T.

        where :math:`h_{Sun}` is the true solar hour angle at the observer's
        longitude and :math:`Eq.T.` is the equation of time.
        """

        if type(self._equation_of_time) == type(None):
            self.equation_of_time(self._julian)

        HRA = self.true_hour_angle(longitude) - self._equation_of_time

        if kwargs:
            HRA = self._interp(HRA, **kwargs)

        return HRA

    def true_solar_time(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return the true solar time (TTs).

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing true solar time in decimal hours.
        
        See Also
        --------
        mean_solar_time

        Notes
        -----
        The true solar time, :math:`TTs`, is given by

        .. math:: TTs = h_{Sun} + 12^h

        where :math:`h_{Sun}` is the true solar hour angle of the at the
        observer's longitude.
        """

        hour_angle = self.true_hour_angle(longitude)
        TST = (hour_angle / 15 + 12) % 24

        if kwargs:
            TST = self._interp(TST, **kwargs)

        return TST

    def mean_solar_time(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return the mean solar time (MTs, same as LMT).

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing mean solar time in decimal hours.
        
        See Also
        --------
        true_solar_time

        Notes
        -----
        The mean solar time, :math:`MTs`, is given by

        .. math:: MTs = h_{mSun} + 12^h

        where :math:`h_{mSun}` is the mean solar hour angle of the at the
        observer's longitude.
        """

        hour_angle = self.mean_hour_angle(longitude)
        MST = (hour_angle / 15 + 12) % 24

        if kwargs:
            MST = self._interp(MST, **kwargs)

        return MST
    
    def sun_right_ascension(self) -> np.array:
        """Return the right ascension of the mean sun position.

        Returns
        -------
        np.array
            Array of shape (n,) containing the right ascension of the mean sun
            in decimal degrees.
        """

        d_U = np.abs(self._julian - 2451545)
        t_U = d_U / 36525
        alpha = (67310.54841 + 8640184.812866 * t_U + 0.093104 * t_U ** 2 - 0.0000062 * t_U ** 3) / 3600

        self._sun_right_ascension = alpha

        return alpha

    def UT1(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return the universal time (same as GMT).

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing universal time in decimal hours.

        Notes
        -----
        The universal time, :math:`UTs` is given by

        .. math:: UTs = h_{mSun}(Gr) + 12^h

        where :math:`h_{mSun}(Gr)` is the mean solar hour angle at
        `longitude=0`.
        """

        if type(self._UT1) != type(None):
            if kwargs:
                return self._interp(self._UT1, **kwargs)
            else:
                return self._UT1

        hour_angle = self.mean_hour_angle(0)
        UT1 = (hour_angle / 15 + 12) % 24
        self._UT1 = UT1

        if kwargs:
            UT1 = self._interp(UT1, **kwargs)

        return UT1
    
    def julian(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return julian times.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing the julian times in decimal days.
        """

        julian = self._julian

        if kwargs:
            julian = self.interp(julian, **kwargs)

        return julian

    def datetime_UTC(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return `datetime.datetime` UTC strings.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing `datetime.datetime` UTC strings.
        """

        jul_data = self._julian

        if kwargs:
            jul_data = self._interp(jul_data, **kwargs)

        datetime = np.zeros((jul_data.size,))
        for i, time in enumerate(jul_data):
            datetime[i] = julian.from_jd(time)

        return datetime
    
    def GMST(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return Greenwhich Mean Sidereal Time.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing Greenwhich Mean Sideral Time in
            decimal hours.
        
        Notes
        -----
        The Greenwhich Mean Sidereal Time can be calculated from the following
        equation:

        .. math:: GMST = \left(280.46061837 + 360.98564736629(j - 2451545) + 0.000387933T^2 - \\frac{T^3}{38710000}\\right) % 360 / 15

        where

        .. math:: T = \\frac{j - 2451545}{36525}

        and :math:`j` is the current Julian day.
        """

        if type(self._GMST) != type(None):
            if kwargs:
                return self._interp(self._GMST, **kwargs)
            else:
                return self._GMST

        julData = self._julian
        T = (julData - 2451545) / 36525
        GMST = (280.46061837 + 360.98564736629 * (julData - 2451545) + 0.000387933 * T **2 - T ** 3 / 38710000) % 360 / 15

        if kwargs:
            GMST = self._interp(GMST, **kwargs)
        
        return GMST
    
    def LMST(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return Local Mean Sidereal Time.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing Local Mean Sideral Time in decimal
            hours.
        
        Notes
        -----
        The Local Mean Sidereal Time can be calculated from the mean solar
        hour angle at the observer's longitude, :math:`h_{mSun}`, and the
        right ascension of the mean Sun position, :math:`\\alpha_{mSun}`, as
        follows:

        .. math:: LMST = h_{mSun} + \\alpha_{mSun}
        """

        if type(self._sun_right_ascension) == type(None):
            self.sun_right_ascension()

        hour_angle = self._mean_hour_angle(longitude) / 15
        alpha = self._sun_right_ascension

        LMST = hour_angle + alpha

        if kwargs:
            LMST = self._interp(LMST, **kwargs)
        
        return LMST
    
    def GAST(self, **kwargs: Dict[str, Any]) -> np.array:
        """Return Greenwhich Apparent Sidereal Time.

        Parameters
        ----------
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing Greenwhich Apparent Sideral Time in
            decimal hours.
        """

        if type(self._GMST) == type(None):
            self.GMST()

        GMST = self._GMST
        EoE = self.equation_of_equinoxes(self._julian)

        GAST = GMST + EoE / 3600

        if kwargs:
            GAST = self.interp(GAST, **kwargs)

        return GAST

    def LAST(self, longitude: np.array, **kwargs: Dict[str, Any]) -> np.array:
        """Return Local Apparent Sidereal Time.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing the actual astronomical longitude
            in decimal degrees.
        kwargs : Dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing Local Apparent Sideral Time in
            decimal hours.
        """

        if type(self._UT1) == type(None):
            self.UT1()

        if type(self._sun_right_ascension) == type(None):
            self.sun_right_ascension()
        
        if type(self._equation_of_equinoxes) == type(None):
            self.equation_of_equinoxes(self._julian)

        UT1 = self._UT1
        alpha_sun = self._sun_right_ascension
        EoE = self._equation_of_equinoxes

        LAST = UT1 - 12 + alpha_sun + EoE + longitude

        if kwargs:
            LAST = self.interp(LAST, **kwargs)

        return LAST
