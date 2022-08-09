

from celest.units.display import _to_string, get_dimension_string


class BaseUnit:
    """BaseUnit()

    The `BaseUnit` class defines the operations required for handling units.
    """

    def __str__(self) -> str:
        return _to_string(self)

    @property
    def dimension(self) -> str:
        """
        Return the dimension string of the unit.
        """
        return get_dimension_string(self)

    @property
    def scale(self) -> float:
        """
        Return the scale of the unit.
        """
        return 1.0

    @property
    def bases(self) -> list:
        """
        Return the list of bases of the unit.
        """
        return [self]

    @property
    def powers(self) -> list:
        """
        Return the list of powers of the unit.
        """
        return [1]

    def __pow__(self, power):
        if not isinstance(power, (int, float)):
            raise TypeError("Units must be raised to scalar powers.")

        return CompoundUnit(1.0, [self], [power])

    def __truediv__(self, other):
        if isinstance(other, BaseUnit):
            return CompoundUnit(1.0, [self, other], [1, -1])
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if other == 1:
            return CompoundUnit(1.0, [self], [-1])
        elif isinstance(other, BaseUnit):
            return CompoundUnit(1.0, [other, self], [1, -1])
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, BaseUnit):
            if other.is_unity():
                return self
            elif self.is_unity():
                return other
            return CompoundUnit(1.0, [self, other], [1, 1])
        elif isinstance(other, (float, int)):
            return CompoundUnit(other, [self], [1])

    def __rmul__(self, other):
        if isinstance(other, BaseUnit):
            if other.is_unity():
                return self
            elif self.is_unity():
                return other
            return CompoundUnit(1.0, [other, self], [1, 1])
        elif isinstance(other, (float, int)):
            return CompoundUnit(other, [self], [1])

    def is_unity(self) -> bool:
        """
        Return `True` if the unit is unity.
        """
        return False


class NamedUnit(BaseUnit):
    """NamedUnit(short_name, long_name, namespace=None)

    The `NamedUnit` class is a `BaseUnit` with a name and associated with a
    namespace.

    Parameters
    ----------
    short_name : str
        The short name of the unit.

        The short name of the unit is used as its namespace identifier.
    long_name : str
        The long name of the unit.
    namespace : dict, optional
        The namespace to insert the unit into.
    """

    def __init__(self, short_name: str, long_name: str, namespace: dict=None) -> None:

        super().__init__()

        self._short_name = short_name
        self._long_name = long_name

        self._insert_into_namespace(namespace)

    def __repr__(self) -> str:
        return "NamedUnit(\"" + _to_string(self) + "\")"

    @property
    def short_name(self) -> str:
        """
        Return the short name of the unit.
        """
        return self._short_name

    @property
    def long_name(self) -> str:
        """
        Return the long name of the unit.
        """
        return self._long_name

    def _insert_into_namespace(self, namespace):

        if namespace is not None and self._short_name not in namespace:
            namespace[self._short_name] = self


class Unit(NamedUnit, BaseUnit):
    """Unit(short_name, long_name, namespace=None, base_units=None)

    The `Unit` class is the main unit object that is associated with a
    namespace and has a name.

    Parameters
    ----------
    short_name : str
        The short name of the unit.

        The short name of the unit is used as its namespace identifier.
    long_name : str
        The long name of the unit.
    namespace : dict, optional
        The namespace to insert the unit into.
    base_units : BaseUnit, optional
        The base units of the unit. That is, how the new unit relates to a base
        SI unit.

        For example, a kilometer unit has the base units of `1000 * u.m`.
    """

    def __init__(self, short_name: str, long_name: str, namespace: dict=None,
                 base_units: BaseUnit=None) -> None:

        super().__init__(short_name, long_name, namespace)

        if base_units is not None:
            self._scale = base_units.scale
            self._bases = base_units.bases
            self._powers = base_units.powers
        else:
            self._scale = 1.0
            self._bases = [self]
            self._powers = [1]

    def __repr__(self):
        return "Unit(\"" + _to_string(self) + "\")"

    @property
    def scale(self) -> float:
        """
        Return the scale of the unit.
        """
        return self._scale

    @property
    def bases(self) -> list:
        """
        Return the list of bases of the unit.
        """
        return self._bases

    @property
    def powers(self) -> list:
        """
        Return the list of powers of the unit.
        """
        return self._powers


class CompoundUnit(BaseUnit):
    """CompoundUnit(scale, bases, powers)

    The `CompoundUnit` class handles more complicated units that are
    combinations of multiple `Unit` objects.

    Parameters
    ----------
    scale : float
        The scale of the unit.
    bases : list
        The list of bases of the unit.
    powers : list
        The list of powers of the unit.
    """

    def __init__(self, scale: float, bases: list, powers: list) -> None:

        self._is_unity = False

        base = bases[0]
        power = powers[0]
        if not isinstance(base, CompoundUnit) and len(bases) == 1 and power >= 0:
            if power == 0:
                self._scale = scale * base.scale
                self._bases = []
                self._powers = []
            else:
                self._scale = scale * base.scale ** power
                self._bases = bases
                self._powers = powers
        else:
            self._scale = scale
            self._bases = bases
            self._powers = powers
            self._expand_and_collect()

    def __repr__(self) -> str:
        return "CompoundUnit(\"" + _to_string(self) + "\")"

    @property
    def scale(self) -> float:
        return self._scale

    @property
    def bases(self) -> list:
        return self._bases

    @property
    def powers(self) -> list:
        return self._powers

    def _expand_and_collect(self) -> None:
        """
        Remove all instances of composite units in `_bases` and `_powers`.
        """

        bases = []
        powers = []

        def handle_repeated_units(repeated_base, repeated_power):
            index = bases.index(repeated_base)
            new_power = powers[index] + repeated_power
            if not new_power:
                bases.pop(index)
                powers.pop(index)
            else:
                powers[index] = new_power

        for base, power in zip(self._bases, self._powers):
            if isinstance(base, CompoundUnit):
                for b, p in zip(base.bases, base.powers):
                    if b not in bases:
                        if b not in bases:
                            bases.append(b)
                            powers.append(p * power)
                        else:
                            handle_repeated_units(b, p * power)
                    else:
                        handle_repeated_units(b, p)
            elif base not in bases:
                bases.append(base)
                powers.append(power)
            else:
                handle_repeated_units(base, power)

        self._bases = bases
        self._powers = powers

        self._is_unity = True if not len(self._bases) else False

    def is_unity(self) -> bool:
        return self._is_unity
