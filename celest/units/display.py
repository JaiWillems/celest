

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


UNIT_TO_DIMENSION_DICTIONARY = {
    "m": "L",
    "mm": "L",
    "cm": "L",
    "km": "L",
    "inch": "L",
    "ft": "L",
    "yd": "L",
    "mi": "L",
    "s": "T",
    "min": "T",
    "hr": "T",
    "dy": "T",
    "jd2000": "D",
    "deg": "A",
    "rad": "A",
    "arcsec": "A",
    "arcmin": "A",
    "hourangle": "A"
}


def get_dimension_string(unit):
    numerator_strings = []
    denominator_strings = []

    dimensions, powers = get_dimensions_from_unit(unit)

    for dimension, power in zip(dimensions, powers):
        unit_string = dimension + ("" if abs(power) == 1 else str(abs(power)))
        if power > 0:
            numerator_strings.append(unit_string)
        else:
            denominator_strings.append(unit_string)

    numerator_string = " ".join(numerator_strings) if len(numerator_strings) else "1"
    denominator_string = " ".join(denominator_strings)

    if not len(denominator_strings):
        return numerator_string
    else:
        return " ".join([numerator_string, "/", denominator_string])


def get_dimensions_from_unit(unit):
    dimensions = []
    dimension_powers = []

    for base, power in zip(unit.bases, unit.powers):
        dimension = UNIT_TO_DIMENSION_DICTIONARY[base.short_name]
        if dimension in dimensions:
            dimension_index = dimensions.index(dimension)
            new_power = dimension_powers[dimension_index] + power
            if new_power == 0:
                dimensions.pop(dimension_index)
                new_power.pop(dimension_index)
            else:
                dimension_powers[dimension_index] = new_power
        else:
            dimensions.append(dimension)
            dimension_powers.append(power)

    return dimensions, dimension_powers
