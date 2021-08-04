# Celest 0.2.0

[![PyPI version](https://badge.fury.io/py/Celest.svg)](https://badge.fury.io/py/Celest) ![PyPI - Downloads](https://img.shields.io/pypi/dm/Celest) ![PyPI - License](https://img.shields.io/pypi/l/celest) [![Documentation Status](https://readthedocs.org/projects/celest/badge/?version=latest)](https://celest.readthedocs.io/en/latest/?badge=latest) ![PyPI - Format](https://img.shields.io/pypi/format/Celest)

The Celest library is designed to provide a simple interface for satellite dynamics and pass analysis calculations.
* **Documentation:** https://celest.readthedocs.io/en/latest/
* **Source Code:** https://github.com/JaiWillems/Celest
* **PyPI:** https://pypi.org/project/Celest/
* **Bug Report or Feature Request:** https://github.com/JaiWillems/Celest/issues

Celest provides:
* Fast orbital conversions between ECI, ECEF, and Horizontal coordinate systems.
* Satellite pass anaysis for encounter planning.
* Supplementary satellite calculations.

## Installation
Celest can be installed from PyPI using the following command:
```terminal
pip install Celest
```

## Release 0.2.2 Features
The Celest 0.2.2 release will include:
* Test driven development,
* Analytical pass analysis method,
* More accurate conversions by accounting for precession, nutation, annd polar motion,
* GCRS and ITRS coordinates and their conversions (replacing ECI and ECEF conversion),
* Epoch selection for obliquity calculations,
* Equinox and solstice time calculations, and
* Moon phase, apogee and perigee calculations.

## Long Term Release Features
Future Celest releases aim to include:
* Encounter optimization and scheduling algorithm for the `Encounter` class,
* Improved runtime using multiprocessing and by incorporating a c-backend, and
* Satellite coordinate conversions for the Equatorial, Ecliptic, Galactic, and Supergalactic coordinate systems with simplicity and efficiency in mind.
