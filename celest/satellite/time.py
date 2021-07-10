"""Time representations."""


from celest.core.decorators import set_module
from celest.satellite import Interpolation
from celest.astronomy import Sun
import numpy as np
import julian


@set_module('celest.satellite')
class Time(object):
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
    equation_of_time(**kwargs)
        Calculate the equation of time.
    true_hour_angle(longitude, **kwargs)
        Calculate the true hour angle.
    mean_hour_angle(longitude, **kwargs)
        Calculate the mean hour angle.
    true_solar_time(longitude, **kwargs)
        Calculate the true solar time.
    mean_solar_time(longitude, **kwargs)
        Calculate the mean solar time.
    UT1(**kwargs)
        Calculate the universal time.
    datetime_UTC(**kwargs)
        Convert the julian times to datetime strings in UTC.
    LMST(longitude, **kwargs)
        Calculate the local mean sidereal time.
    GMST(longitude, **kwargs)
        Calculate the Greenwhich mean sidereal time.
    """

    def __init__(self, julian: np.array, offset: float=0, factor: int=0) -> None:
        """Initialize attributes."""

        self._equation_of_time = None
        self._UT1 = None
        self._GMST = None
        self._GAST = None

        self._interp = Interpolation()

        if factor > 0:
            julian = self._interp(data=julian, factor=factor)

        self.julian = julian + offset
        self.length = self.julian.size

    def equation_of_time(self, **kwargs) -> np.array:
        """Calculate the equation of time.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing the equation of time in decimal
            minutes.

        Notes
        -----
        This method uses an empirical equation to approximate the equation of
        time within an accuracy of 0.5 minutes.

        .. math:: EoT = 9.87\sin(2B) - 7.53\cos(B) - 1.5\sin(B)

        where :math:`B = \\frac{360}{365}(d - 81)` and :math:`b` is the number
        of days since the start of the year.
        """

        if type(self._equation_of_time) != type(None):
            if kwargs:
                return self._interp(self._equation_of_time)
            else:
                return self._equation_of_time

        day_num = np.zeros((self.length,))
        for i, time in enumerate(self.julian):
            day_num[i] = float(julian.from_jd(time).strftime("%j"))

        B = (day_num - 81) * 360 / 365
        B_r, B2_r = np.radians(B), np.radians(2 * B)
        EoT = 9.87 * np.sin(B2_r) - 7.53 * np.cos(B_r) - 1.5 * np.sin(B_r)

        self._equation_of_time = EoT

        if kwargs:
            EoT = self._interp(EoT, **kwargs)

        return EoT

    def true_hour_angle(self, longitude: np.array, **kwargs) -> np.array:
        """Calculate the solar hour angle using the true sun position.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : dict, optional
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

    def mean_hour_angle(self, longitude: np.array, **kwargs) -> np.array:
        """Calculate the hour angle of the mean solar position.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : dict, optional
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
            self.equation_of_time()

        HRA = self.true_hour_angle(longitude) - self._equation_of_time

        if kwargs:
            HRA = self._interp(HRA, **kwargs)

        return HRA

    def true_solar_time(self, longitude: np.array, **kwargs) -> np.array:
        """Get the true solar time (TTs).

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : dict, optional
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

    def mean_solar_time(self, longitude: np.array, **kwargs) -> np.array:
        """Get the mean solar time (MTs, same as LMT).

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : dict, optional
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
    
    def _alpha_mean_sun(self) -> np.array:
        """Get the right ascension of the mean sun position.

        Returns
        -------
        np.array
            Array of shape (n,) containing the right ascension of the mean sun
            in decimal degrees.
        """

        d_U = np.abs(self.julian - 2451545)
        t_U = d_U / 36525
        alpha = 67310.54841 + 8640184.812866 * t_U + 0.093104 * t_U ** 2 - 0.0000062 * t_U ** 3

        return alpha / 3600

    def UT0(self, longitude: np.array, **kwargs) -> np.array:
        pass

    def UT1(self, **kwargs) -> np.array:
        """Get the universal time (same as GMT).

        Parameters
        ----------
        kwargs : dict, optional
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

    def datetime_UTC(self, **kwargs) -> np.array:
        """Get `datetime.datetime` UTC strings

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing `datetime.datetime` UTC strings.
        """

        if kwargs:
            jul_data = self._interp(self.julian, **kwargs)
        else:
            jul_data = self.julian

        datetime = np.zeros((jul_data.size,))
        for i, time in enumerate(jul_data):
            datetime[i] = julian.from_jd(time)

        return datetime

    def LMST(self, longitude: np.array, **kwargs) -> np.array:
        """Get local mean sidereal time.

        Parameters
        ----------
        longitude : np.array
            Array of shape (n,) containing longitude in decimal degrees.
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing local mean sideral time in decimal
            hours.
        
        Notes
        -----
        The local mean sidereal time can be calculated from the mean solar
        hour angle at the observer's longitude, :math:`h_{mSun}`, and the
        right ascension of the mean Sun position, :math:`\\alpha_{mSun}`, as
        follows:

        .. math:: LMST = h_{mSun} + \\alpha_{mSun}
        """

        hour_angle = self._mean_hour_angle(longitude) / 15
        alpha = self._alpha_mean_sun()

        LMST = hour_angle + alpha

        if kwargs:
            LMST = self._interp(LMST, **kwargs)
        
        return LMST

    def LAST(self, longitude: np.array, **kwargs) -> np.array:
        pass

    def GMST(self, **kwargs) -> np.array:
        """Get Greenwhich mean sidereal time.

        Parameters
        ----------
        kwargs : dict, optional
            Optional keyword arguments passed into the `Interpolation()`
            callable.

        Returns
        -------
        np.array
            Array of shape (n,) containing greenwhich mean sideral time in
            decimal hours.
        
        Notes
        -----
        The Greenwhich mean sidereal time can be calculated from the mean
        solar hour angle at `longitude=0`, :math:`h_{mSun}(Gr)`, and the right
        ascension of the mean Sun position, :math:`\\alpha_{mSun}`, as follows:

        .. math:: GMST = h_{mSun}(Gr) + \\alpha_{mSun}
        """

        if type(self._GMST) != type(None):
            if kwargs:
                return self._interp(self._GMST, **kwargs)
            else:
                return self._GMST
        
        hour_angle = self._mean_hour_angle(0) / 15
        alpha = self._alpha_mean_sun()

        GMST = hour_angle + alpha
        self._GMST = GMST

        if kwargs:
            GMST = self._interp(GMST, **kwargs)
        
        return GMST

    def GAST(self, **kwargs) -> np.array:
        pass
