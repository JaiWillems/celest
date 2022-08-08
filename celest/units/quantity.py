

from celest.units.core import CompoundUnit
from copy import deepcopy


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
        Return a new Quantity object with the new unit.
    _convert_to(new_unit)
        Convert Quantity data to a different unit.


    Examples
    --------
    >>> q = Quantity(5, u.m)
    """

    def __init__(self, data, unit):

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
    def data(self):
        return self._data

    @property
    def unit(self):
        return self._unit

    @property
    def dimension(self):
        return self._unit.dimension

    def to(self, new_unit):
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

        return self._convert_to(new_unit).data

    def _convert_to(self, new_unit):
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

        >>> quantity_in_kilometers = quantity_in_meters._convert_to(u.km)
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
