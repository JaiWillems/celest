

from celest.units.utils import setup_unit


namespace = globals()


# Si length measures.
setup_unit("m", "meter", namespace)
setup_unit("mm", "millimeter", namespace, 0.001 * m)
setup_unit("cm", "centimeter", namespace, 0.01 * m)
setup_unit("km", "kilometer", namespace, 1000 * m)


# Imperial length measures.
setup_unit("inch", "inch", namespace, 0.0254 * m)
setup_unit("ft", "feet", namespace, 0.305 * m)
setup_unit("yd", "yard", namespace, 0.914 * m)
setup_unit("mi", "mile", namespace, 1609.344 * m)


# Time measures.
setup_unit("s", "second", namespace)
setup_unit("min", "minute", namespace, 60 * s)
setup_unit("hr", "hour", namespace, 3600 * s)
setup_unit("dy", "day", namespace, 86400 * s)


# Date measures.
setup_unit("jd2000", "Julian day 2000", namespace)


# Angular measures.
setup_unit("deg", "degree", namespace)
setup_unit("rad", "radian", namespace, 57.29577951308232 * deg)
setup_unit("arcsec", "seconds of arc", namespace, (1 / 3600) * deg)
setup_unit("arcmin", "minutes of arc", namespace, (1 / 60) * deg)
setup_unit("hourangle", "hour angle", namespace, 15 * deg)
