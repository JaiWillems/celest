# <a href="https://celest.readthedocs.io/en/latest/"><img alt="Celest" src="/branding/logo/transparent_background_large.png" height="150"></a>

[![PyPI version](https://badge.fury.io/py/celest.svg)](https://badge.fury.io/py/celest) [![license](https://img.shields.io/pypi/l/celest)](https://img.shields.io/pypi/l/celest) [![Downloads](https://static.pepy.tech/personalized-badge/celest?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=downloads)](https://pepy.tech/project/celest) [![Documentation Status](https://readthedocs.org/projects/celest/badge/?version=latest)](https://celest.readthedocs.io/en/latest/?badge=latest) [![python](https://img.shields.io/pypi/pyversions/celest)](https://img.shields.io/pypi/pyversions/celest) [![format](https://img.shields.io/pypi/wheel/celest)](https://img.shields.io/pypi/wheel/celest)

Celest is a satellite dynamics and mission planning library designed for the University of Toronto Aerospace Team which is applying the program to a hyperspectral imaging CubeSat mission. The library aims to provide tools to plan satellite-ground interactions such as ground target tracking, and calculating imaging and transmission opportunities.

* **Documentation:** https://celest.readthedocs.io/en/latest/
* **Source Code:** https://github.com/JaiWillems/celest
* **PyPI:** https://pypi.org/project/celest/
* **Bug Report or Feature Request:** https://github.com/JaiWillems/Celest/issues

## Installation

Celest can be installed from PyPI using the command line:

```terminal
pip install celest
```

## What can the program do?

The program requires an input data set containing satellite positions in the gcrs, itrs, or geographical frame with corresponding Julian times. The program can then calculate various position and time representations, satellite to ground encounters, and transformations necessary for ground location tracking.

## What can be expected of the project?

In future releases, a high order c-implemented orbit propagator will be incorporated into the library to allow for a complete mission planning workflow from inputting orbital parameters to calculating satellite-to-ground encounter opportunities.

We also intend to incorporate a window scheduling algorithm to create an optimal satellite itinerary for satellite mission planning.

Other features under consideration include a terminal interface, higher fidelity conversions, window/data visuals, data loaders, and more robust ground tracking by incorporating forward motion compensation.
