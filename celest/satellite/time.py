

from celest.core.decorators import set_module
from celest.satellite import Interpolation
from celest.astronomy import Sun
from typing import Union
import numpy as np
import julian


@set_module('celest.satellite')
class Time(object):

    def __init__(self, julian: Union[float, np.array], offset: float=0, factor:
                 int=0) -> None:

        self._equation_of_time = None
        self._UT1 = None
        self._GMST = None
        self._GAST = None

        self.interp = Interpolation()

        if factor > 0:
            julian = self.interp(data=julian, factor=factor)

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
                return self.interp(self._equation_of_time)
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
            EoT = self.interp(EoT, **kwargs)

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

        sun_dot_nb = np.diag(np.matmul(sun_position, norm_b.T))
        denom = np.linalg.norm(norm_b, axis=1) ** 2
        scale_arr = np.repeat((sun_dot_nb / denom).reshape((-1, 1)), 3, axis=1)
        proj_v4_on_nb = scale_arr * norm_b

        test = np.diag(np.matmul(proj_v4_on_nb, norm_b.T))

        pos_ind = np.where(test >= 0)[0]

        HRA[pos_ind] = -HRA[pos_ind]

        if kwargs:
            HRA = self.interp(HRA, **kwargs)

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
        """

        if type(self._equation_of_time) == type(None):
            self.equation_of_time()

        HRA = self.true_hour_angle(longitude) - self._equation_of_time

        if kwargs:
            HRA = self.interp(HRA, **kwargs)

        return HRA

    def true_solar_time(self, longitude: np.array, **kwargs) -> np.array:
        """Get the true solar time.

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
        """

        hour_angle = self.true_hour_angle(longitude)
        TST = (hour_angle / 15 + 12) % 24

        if kwargs:
            TST = self.interp(TST, **kwargs)

        return TST

    def mean_solar_time(self, longitude: np.array, **kwargs) -> np.array:
        """Get the mean solar time.

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
        """

        hour_angle = self.mean_hour_angle(longitude)
        MST = (hour_angle / 15 + 12) % 24

        if kwargs:
            MST = self.interp(MST, **kwargs)

        return MST

    def UT0(self, longitude: np.array, **kwargs) -> np.array:
        pass

    def UT1(self, **kwargs) -> np.array:
        pass

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
            jul_data = self.interp(self.julian, **kwargs)
        else:
            jul_data = self.julian

        datetime = np.zeros((jul_data.size,))
        for i, time in enumerate(jul_data):
            datetime[i] = julian.from_jd(time)

        return datetime

    def LMST(self, longitude: np.array, **kwargs) -> np.array:
        pass

    def LAST(self, longitude: np.array, **kwargs) -> np.array:
        pass

    def GMST(self, **kwargs) -> np.array:
        pass

    def GAST(self, **kwargs) -> np.array:
        pass