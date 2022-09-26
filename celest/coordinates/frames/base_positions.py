

from celest.units.core import Unit
from celest.units.quantity import Quantity
import numpy as np


class Position2d:
    """Position2d(x, x_unit, y, y_unit, time, time_unit)

    2D position data container.

    Parameters
    ----------
    x, y, time : np.ndarray
        The spatial and temporal coordinate dimensions.
    x_unit, y_unit, time_unit : Unit
        Units of the spatial and temporal coordinate dimensions.
    """

    def __init__(self, x: np.ndarray, x_unit: Unit, y: np.ndarray, y_unit:
                 Unit, time: np.ndarray, time_unit: Unit) -> None:

        if time.ndim != 1 or x.ndim != 1 or y.ndim != 1:
            raise ValueError("Input arrays should be one dimensional.")
        if time.size != x.size != y.size:
            raise ValueError("Input arrays should have the same length.")

        self._x = Quantity(x, x_unit)
        self._y = Quantity(y, y_unit)
        self._time = Quantity(time, time_unit)

    def _get_x(self) -> Quantity:
        return self._x

    def _get_y(self) -> Quantity:
        return self._y

    @property
    def time(self) -> Quantity:
        return self._time

    def to_numpy(self, unit: Unit) -> np.ndarray:
        """Convert coordinate data to NumPy array.

        Parameters
        ----------
        unit : Unit
            The unit of the output array.

        Returns
        -------
        np.ndarray
            The spatial coordinate data in the specified unit.
        """

        return np.array([
            self._x.to(unit).data,
            self._y.to(unit).data
        ]).T


class Position3d:
    """Position2d(x, x_unit, y, y_unit, z, z_unit, time, time_unit)

    3D position data container.

    Parameters
    ----------
    x, y, z, time : np.ndarray
        The spatial and temporal coordinate dimensions.
    x_unit, y_unit, z_unit time_unit : Unit
        Units of the spatial and temporal coordinate dimensions.
    """

    def __init__(self, x: np.ndarray, x_unit: Unit, y: np.ndarray, y_unit:
                 Unit, z: np.ndarray, z_unit: Unit, time: np.ndarray,
                 time_unit: Unit) -> None:

        if time.ndim != 1 or x.ndim != 1 or y.ndim != 1 or z.ndim != 1:
            raise ValueError("Input arrays should be one dimensional.")
        if time.size != x.size != y.size != z.size:
            raise ValueError("Input arrays should have the same length.")

        self._x = Quantity(x, x_unit)
        self._y = Quantity(y, y_unit)
        self._z = Quantity(z, z_unit)
        self._time = Quantity(time, time_unit)

    def _get_x(self) -> Quantity:
        return self._x

    def _get_y(self) -> Quantity:
        return self._y

    def _get_z(self) -> Quantity:
        return self._z

    @property
    def time(self) -> Quantity:
        return self._time

    def to_numpy(self, unit: Unit) -> np.ndarray:
        """Convert coordinate data to NumPy array.

        Parameters
        ----------
        unit : Unit
            The unit of the output array.

        Returns
        -------
        np.ndarray
            The spatial coordinate data in the specified unit.
        """

        return np.array([
            self._x.to(unit).data,
            self._y.to(unit).data,
            self._z.to(unit).data
        ]).T
