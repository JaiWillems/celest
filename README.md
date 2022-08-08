# <a href="https://celest.readthedocs.io/en/latest/"><img alt="Celest" src="/branding/logo/celest_logo_transparent_wide.png" height="150"></a>

[![PyPI version](https://badge.fury.io/py/celest.svg)](https://pypi.org/project/celest/) [![license](https://img.shields.io/pypi/l/celest)](https://github.com/JaiWillems/celest/blob/main/LICENSE.txt) [![Downloads](https://static.pepy.tech/personalized-badge/celest?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=downloads)](https://pepy.tech/project/celest) [![Documentation Status](https://readthedocs.org/projects/celest/badge/?version=latest)](https://celest.readthedocs.io/en/latest/?badge=latest) [![python](https://img.shields.io/pypi/pyversions/celest)](https://img.shields.io/pypi/pyversions/celest) [![format](https://img.shields.io/pypi/wheel/celest)](https://img.shields.io/pypi/wheel/celest)

Celest is a satellite mission planning software designed for remote sensing CubeSat mission proifles. The library aims to provide the necessary tools for celestial orbital conversions and satellite-to-ground interactions such as determining and scheduling imaging/transmission opportunities.

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

Other features under consideration include a terminal interface, higher fidelity conversions, window/data visuals, data loaders, and more robust ground tracking by incorporating forward motion compensation.
