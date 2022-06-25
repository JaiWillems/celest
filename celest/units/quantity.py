

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

    def to(self, unit):
        """Return data with the new unit.

        Parameters
        ----------
        unit : Unit, CompoundUnit

        Returns
        -------
        Any
            Data in the new unit.

        Examples
        --------
        >>> quantity_in_kilometers = quantity_in_meters.to(u.km)
        """
        if self._unit.dimension != unit.dimension:
            raise ValueError("Unit dimensionality is mismatched.")
        scale = self._get_unit_scale(self.unit) / self._get_unit_scale(unit)
        return self._data * scale

    def _get_unit_scale(self, unit):
        if isinstance(unit, CompoundUnit):
            scale = 1.0
            for base, power in zip(unit.bases, unit.powers):
                scale *= base.scale ** power
            return scale
        else:
            return unit.scale

    def get_unit(self):
        return self._unit

    def convert_to(self, unit):
        """Convert data to different units.

        Parameters
        ----------
        unit : Unit, CompoundUnit
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

        if self._unit.dimension != unit.dimension:
            raise ValueError("Unit dimensionality is mismatched.")

        self._data *= self._get_unit_scale(self.unit) / self._get_unit_scale(unit)
        self._unit = unit
