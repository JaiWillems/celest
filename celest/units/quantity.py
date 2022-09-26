

from celest.units.core import CompoundUnit, Unit
from copy import deepcopy
from typing import Any, Union


class Quantity:
    """Quantity(data, unit)

    Data with a unit.

    The Quantity class holds both data and a unit together to allow for easy
    unit conversions.

    Parameters
    ----------
    data : Any
        The data to be stored in the Quantity object.

        To allow for unit conversions, the data object must have the
        multiplication operations defined.
    unit : Unit, CompoundUnit
        The unit associated with the data.

    Attributes
    ----------
    data : Any
        The data stored in the Quantity.
    unit : Unit, CompoundUnit
        The unit associated with the data.

    Methods
    -------
    to(new_unit)
        Return the Quantity data in the new unit.

    Examples
    --------
    >>> q = Quantity(5, u.m)

    The data data, unit, and dimension string are available as attributes.

    >>> q1.data
    5
    >>> q1.unit
    Unit("m")
    >>> q1.dimension
    'L'

    Quantities of the same dimension can be added or subtracted.

    >>> q2 = Quantity(0.1, u.km)
    >>> q3 = q1 + q2
    >>> print(q3)
    105.0 m

    >>> q4 = q1 - q2
    >>> print(q4)
    -95.0 m

    Quantities can be multiplied or divided.

    >>> q5 = Quantity(2, u.s)
    >>> q6 = q1 * q5
    >>> print(q6)
    10 m s

    >>> q7 = q1 / q5
    >>> print(q7)
    2.5 m / s
    """

    def __init__(self, data: Any, unit: Union[Unit, CompoundUnit]) -> None:
        """Data with a unit.

        The Quantity class holds both data and a unit together to allow for easy
        unit conversions.

        Parameters
        ----------
        data : Any
            The data to be stored in the Quantity object.

            To allow for unit conversions, the data object must have the
            multiplication operations defined.
        unit : Unit, CompoundUnit
            The unit associated with the data.
        """

        self._data = data
        self._unit = unit

    def __str__(self):
        return str(self._data) + " " + str(self._unit)

    def __repr__(self):
        return "Quantity(" + repr(self._data) + ", " + repr(self._unit) + ")"

    def __neg__(self):
        return Quantity(-self._data, self._unit)

    def __add__(self, other):
        if isinstance(other, (float, int)):
            data = self._data + other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            if self._unit.dimension != other._unit.dimension:
                raise ArithmeticError(
                    "Unit dimensions must match for addition/subtraction.")
            data = self._data + other.to(self._unit)
            return Quantity(data, self._unit)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return -self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other, (float, int)):
            data = self.data * other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            if self._unit.dimension == other._unit.dimension:
                data = self._data * other.to(self._unit)
                unit = self._unit ** 2
            else:
                data = self._data * other._data
                unit = self._unit * other._unit
            return Quantity(data, unit)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (float, int)):
            data = self._data / other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            if self._unit.dimension == other._unit.dimension:
                data = self._data / other.to(self._unit)
                unit = self._unit / self._unit
            else:
                data = self._data / other._data
                unit = self._unit / other._unit
            return Quantity(data, unit)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, (float, int)):
            data = other / self._data
            unit = 1 / self._unit
            return Quantity(data, unit)
        elif isinstance(other, Quantity):
            if self._unit.dimension == other._unit.dimension:
                data = other.to(self._unit) / self._data
                unit = self._unit / self._unit
            else:
                data = other._data / self._data
                unit = other._unit / self._unit
            return Quantity(data, unit)
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Quantity):
            if self._unit.dimension != other._unit.dimension:
                raise ArithmeticError(
                    "Unit dimensions must match for comparison.")
            return self._data == other.to(self._unit)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, Quantity):
            if self._unit.dimension != other._unit.dimension:
                raise ArithmeticError(
                    "Unit dimensions must match for comparison.")
            return self._data < other.to(self._unit)
        else:
            return NotImplemented

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __gt__(self, other):
        if isinstance(other, Quantity):
            if self._unit.dimension != other._unit.dimension:
                raise ArithmeticError(
                    "Unit dimensions must match for comparison.")
            return self._data > other.to(self._unit)
        else:
            return NotImplemented

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)

    @property
    def data(self) -> Any:
        """
        Return the data of the quantity.
        """
        return self._data

    @property
    def unit(self) -> Union[Unit, CompoundUnit]:
        """
        Return the unit of the quantity.
        """
        return self._unit

    @property
    def dimension(self) -> str:
        """
        Return the dimension string of the quantity.
        """
        return self._unit.dimension

    def to(self, new_unit: Union[Unit, CompoundUnit]) -> Any:
        """Return the Quantity data in the new unit.

        Parameters
        ----------
        new_unit : Unit, CompoundUnit
            New unit with the same dimension as the current unit.

        Returns
        -------
        Any
            The Quantity data in the new unit.

        Examples
        --------
        >>> quantity_in_meters = Quantity(5, u.m)
        >>> quantity_in_meters.to(u.km)
        0.005
        """

        return self.convert_to(new_unit).data

    def convert_to(self, new_unit: Union[Unit, CompoundUnit]) -> "Quantity":
        """Return a new Quantity object with the new unit.

        Parameters
        ----------
        new_unit : Unit, CompoundUnit
            New unit with the same dimension as the current unit.

        Returns
        -------
        Quantity
            A new Quantity object with the new unit.

        Examples
        --------
        Initialize the `Quantity` object:

        >>> quantity_in_meters = Quantity(5, u.m)

        Convert to kilometers:

        >>> quantity_in_kilometers = quantity_in_meters.convert_to(u.km)
        """

        if self._unit.dimension != new_unit.dimension:
            raise ValueError(
                "New unit has a different dimension then the current unit.")

        def _get_unit_scale(unit):
            if isinstance(unit, CompoundUnit):
                scale = 1.0
                for base, power in zip(unit.bases, unit.powers):
                    scale *= base.scale ** power
                return scale
            else:
                return unit.scale

        quantity_copy = deepcopy(self)

        if repr(quantity_copy._unit) != repr(new_unit):
            quantity_copy._data *= _get_unit_scale(
                quantity_copy.unit) / _get_unit_scale(new_unit)
            quantity_copy._unit = new_unit

        return quantity_copy
