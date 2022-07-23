

from celest.units.core import CompoundUnit
from copy import deepcopy


class Quantity:
    """Quantity(data, unit)

    Data with a unit.

    Parameters
    ----------
    data : Any
    unit : Unit, CompoundUnit
        The unit associated with the data.

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
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError(
                    "Units must match for Quantity addition.")
            data = self._data + other._data
            return Quantity(data, self._unit)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, (float, int)):
            data = self._data + other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError(
                    "Units must match for Quantity addition.")
            data = self._data + other._data
            return Quantity(data, self._unit)
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (float, int)):
            data = self._data - other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError(
                    "Units must match for Quantity subtraction.")
            data = self._data - other._data
            return Quantity(data, self._unit)
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, (float, int)):
            data = other - self._data
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError(
                    "Units must match for Quantity subtraction.")
            data = other._data - self._data
            return Quantity(data, self._unit)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (float, int)):
            data = self.data * other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            data = self._data * other._data
            unit = self._unit * other._unit
            return Quantity(data, unit)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (float, int)):
            data = self._data * other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
            data = self._data * other._data
            unit = self._unit * other._unit
            return Quantity(data, unit)
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, (float, int)):
            data = self._data / other
            return Quantity(data, self._unit)
        elif isinstance(other, Quantity):
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
            data = other._data / self._data
            unit = other._unit / self._unit
            return Quantity(data, unit)
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError("Units must match for comparison.")
            return self._data == other._data
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError("Units must match for comparison.")
            return self._data != other._data
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError("Units must match for comparison.")
            return self._data < other._data
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError("Units must match for comparison.")
            return self._data <= other._data
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError("Units must match for comparison.")
            return self._data > other._data
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, Quantity):
            if repr(self._unit) != repr(other._unit):
                raise ArithmeticError("Units must match for comparison.")
            return self._data >= other._data
        else:
            return NotImplemented

    @property
    def data(self):
        return self._data

    @property
    def unit(self):
        return self._unit

    def to(self, new_unit):
        """Return Quantity with the new unit.

        Parameters
        ----------
        new_unit : Unit, CompoundUnit

        Returns
        -------
        Quantity

        Examples
        --------
        >>> quantity_in_kilometers = quantity_in_meters.to(u.km)
        """

        quantity_in_new_unit = deepcopy(self)
        quantity_in_new_unit.convert_to(new_unit)
        return quantity_in_new_unit

    def get_unit(self):
        return self._unit

    def convert_to(self, new_unit):
        """Convert data to different units.

        Parameters
        ----------
        new_unit : Unit, CompoundUnit
            New unit with the same dimension as the current unit.

        Examples
        --------
        Initialize the `Quantity` object:

        >>> q = Quantity(5, u.m)
        >>> q.data
        5

        Convert to kilometers:

        >>> q.convert_to(u.km)
        >>> q.data
        0.005
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

        if repr(self._unit) == repr(new_unit):
            return self
        else:
            self._data *= _get_unit_scale(self.unit) / _get_unit_scale(new_unit)
            self._unit = new_unit
