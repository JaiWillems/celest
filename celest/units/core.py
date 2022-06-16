

from celest.units.display import _to_string


class BaseUnit:
    """BaseUnit()

    The `BaseUnit` class defines the operations required for handling units.
    """

    def __str__(self):
        return _to_string(self)

    @property
    def scale(self):
        """
        Return scale of the unit.
        """
        return 1.0

    @property
    def bases(self):
        """
        Return the bases of the unit.
        """
        return [self]

    @property
    def powers(self):
        """
        Return the powers of the unit.
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
        if isinstance(other, BaseUnit):
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

    def is_unity(self):
        """
        Return `True` if the unit is unity.
        """
        return False


class NamedUnit(BaseUnit):
    """NamedUnit(short_name, long_name, namespace=None)

    The `NamedUnit` class is a `BaseUnit` with a name and associated with a
    namespace.
    """

    def __init__(self, short_name, long_name, namespace=None):

        super().__init__()

        self._short_name = short_name
        self._long_name = long_name

        self.insert_into_namespace(namespace)

    def __repr__(self):
        return "NamedUnit(\"" + _to_string(self) + "\")"

    @property
    def short_name(self):
        return self._short_name

    @property
    def long_name(self):
        return self._long_name

    def insert_into_namespace(self, namespace):

        if namespace is not None and self._short_name not in namespace:
            namespace[self._short_name] = self


# TODO: Is there ever a need to not pass in a namespace?
class Unit(NamedUnit, BaseUnit):
    """Unit(short_name, long_name, namespace=None, base_units=None)

    The `Unit` class is the main unit object that is associated with a
    namespace and has a name.
    """

    def __init__(self, short_name, long_name, namespace=None, base_units=None):

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
    def scale(self):
        """
        Return scale of the unit.
        """
        return self._scale

    @property
    def bases(self):
        """
        Return the bases of the unit.
        """
        return self._bases

    @property
    def powers(self):
        """
        Return the powers of the unit.
        """
        return self._powers


class CompoundUnit(BaseUnit):
    """CompoundUnit(scale, bases, powers)

    The `CompoundUnit` class handles more complicated units that are
    combinations of multiple `Unit` objects.
    """

    def __init__(self, scale, bases, powers):

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

    def __repr__(self):
        return "CompoundUnit(\"" + _to_string(self) + "\")"

    @property
    def scale(self):
        return self._scale

    @property
    def bases(self):
        return self._bases

    @property
    def powers(self):
        return self._powers

    def _expand_and_collect(self):
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
                        bases.append(b)
                        powers.append(p * power)
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

    def is_unity(self):
        return self._is_unity
