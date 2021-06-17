# Celest

[![PyPI version](https://badge.fury.io/py/Celest.svg)](https://badge.fury.io/py/Celest) ![PyPI - Downloads](https://img.shields.io/pypi/dm/Celest) ![PyPI - License](https://img.shields.io/pypi/l/celest) [![Documentation Status](https://readthedocs.org/projects/celest/badge/?version=latest)](https://celest.readthedocs.io/en/latest/?badge=latest) ![PyPI - Format](https://img.shields.io/pypi/format/Celest)

The Celest library is designed to provide a simple interface for satellite positional representations and encounter planning.
* **Documentation:** https://celest.readthedocs.io/en/latest/
* **Source Code:** https://github.com/JaiWillems/Celest
* **PyPI:** https://pypi.org/project/Celest/
* **Bug Report or Feature Request:** https://github.com/JaiWillems/Celest/issues

Celest provides:
* Fast orbital conversions between ECI, ECEF, and Horizontal coordinate systems.
* Encounter generation and planning.

## Installation
Celest can be installed from PyPI with the following command:
```terminal
pip install Celest
```


## Release 0.2.0 Features
The Celest 0.2.0 release will include:
* Test driven development,
* Improved specific and general encounter statistics,
* Sun-encounter constraint angle,
* Analytical pass analysis method,
* ECEF and Geographical coordinate conversions,
* Sexagesimal formatting option for angular outputs,
* Astronomy module for celestial object localization,
* Interpolation-factor parameter for inputed data,
* Public special interpolation methods for all position representations,
* Method for encounter indices return, and
* Faster ECI and ECEF conversions through new multiplication method.

## Long Term Release Features
Goals for future Celest features include:
* Encounter optimization and scheduling algorithm for the Encounter class,
* Improved runtime using multiprocessing, and
* Satellite coordinate conversions for the Equatorial, Ecliptic, Galactic, and Supergalactic coordinate systems with simplicity and efficiency in mind.
