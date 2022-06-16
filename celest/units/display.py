

from celest.units import core  # Import core to prevent circular import errors.


def _to_string(unit):

    if isinstance(unit, core.CompoundUnit):
        return generate_compound_unit_string(unit)
    elif isinstance(unit, core.NamedUnit):
        return unit.short_name
    else:
        raise TypeError("Unit type has no string representation.")


def generate_compound_unit_string(unit):

    numerator_strings = []
    denominator_strings = []
    for base, power in zip(unit.bases, unit.powers):
        unit_string = base.short_name + (
            "" if abs(power) == 1 else str(abs(power)))
        if power > 0:
            numerator_strings.append(unit_string)
        else:
            denominator_strings.append(unit_string)

    scale_string = "" if unit.scale == 1 else (str(unit.scale) + " ")
    numerator_string = " ".join(numerator_strings) if len(
        numerator_strings) else "1"
    denominator_string = " ".join(denominator_strings)

    if not len(denominator_strings):
        return scale_string + numerator_string
    else:
        return " ".join(
            [scale_string + numerator_string, "/", denominator_string])
