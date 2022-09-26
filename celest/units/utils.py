

from celest.units.core import Unit


def setup_unit(short_name, long_name, namespace, base_units=None):
    Unit(short_name, long_name, namespace, base_units)
