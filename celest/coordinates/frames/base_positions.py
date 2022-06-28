

from celest.units.quantity import Quantity


class Position2d:
    """2D position data container.

    Parameters
    ----------
    x, y : Quantity
        The x and y coordinate dimensions.
    time : Quantity, optional
        The time coordinate dimension.
    """

    def __init__(self, x: Quantity, y: Quantity, time: Quantity=None) -> None:

        if not isinstance(x, Quantity) or not isinstance(y, Quantity):
            raise ValueError("x and y should be Quantity objects.")
        if time is not None and not isinstance(time, Quantity):
            raise ValueError("time should be a Quantity object.")

        self._x = x
        self._y = y
        self._time = time

    @property
    def x(self) -> Quantity:
        return self._x

    @property
    def y(self) -> Quantity:
        return self._y

    @property
    def time(self) -> Quantity:
        return self._time


class Position3d:
    """3D position data container.

        Parameters
        ----------
        x, y, z : Quantity
            The x and y coordinate dimensions.
        time : Quantity, optional
            The time coordinate dimension.
        """

    def __init__(self, x: Quantity, y: Quantity, z: Quantity, time: Quantity=None) -> None:

        if not isinstance(x, Quantity) or not isinstance(y, Quantity) or not isinstance(z, Quantity):
            raise ValueError("x, y, and z should be Quantity objects.")
        if time is not None and not isinstance(time, Quantity):
            raise ValueError("time should be a Quantity object.")

        self._x = x
        self._y = y
        self._z = z
        self._time = time

    @property
    def x(self) -> Quantity:
        return self._x

    @property
    def y(self) -> Quantity:
        return self._y

    @property
    def z(self) -> Quantity:
        return self._z

    @property
    def time(self) -> Quantity:
        return self._time
