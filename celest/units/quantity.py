

from celest.units.core import CompoundUnit


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

    @property
    def data(self):
        return self._data

    @property
    def unit(self):
        return self._unit

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
            raise ValueError("Unit dimensionality is mismatched.")

        def _get_unit_scale(unit):
            if isinstance(unit, CompoundUnit):
                scale = 1.0
                for base, power in zip(unit.bases, unit.powers):
                    scale *= base.scale ** power
                return scale
            else:
                return unit.scale

        self._data *= _get_unit_scale(self.unit) / _get_unit_scale(new_unit)
        self._unit = new_unit
