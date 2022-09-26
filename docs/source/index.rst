.. Celest documentation master file, created by
   sphinx-quickstart on Mon Jun  7 18:56:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: ../../branding/logo/celest_logo_transparent_wide.png
   :width: 400
   :alt: Celest Logo

.. image:: https://badge.fury.io/py/celest.svg
   :target: https://pypi.org/project/celest/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/l/celest
   :target: https://github.com/JaiWillems/celest/blob/main/LICENSE.txt
   :alt: License

.. image:: https://static.pepy.tech/personalized-badge/celest?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=downloads
   :target: https://pepy.tech/project/celest
   :alt: Downloads

.. image:: https://readthedocs.org/projects/celest/badge/?version=latest
   :target: https://celest.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/pyversions/celest
   :target: https://img.shields.io/pypi/pyversions/celest
   :alt: Python

.. image:: https://img.shields.io/pypi/wheel/celest
   :target: https://img.shields.io/pypi/wheel/celest
   :alt: format

Celest Overview
---------------

Celest is a satellite mission planning software designed for a remote sensing
CubeSat mission profile. The library aims to provide the necessary tools for
celestial orbital conversions and satellite-to-ground interactions such as
determining and scheduling imaging/transmission opportunities.

What can be expected of the project?
------------------------------------

In future releases, a high order c-implemented orbit propagator will be
incorporated into the library to allow for a complete mission planning workflow
from inputting orbital parameters to calculating satellite-to-ground encounter
opportunities.

Other features that are being considered include a terminal interface, higher
fidelity conversions, window/data visuals, data loaders, and more
robust ground tracking by incorporating forward motion compensation.

Index
-----

.. toctree::
   :maxdepth: 2
   :caption: Guide

   guide/quickstart
   guide/tutorials
   guide/contributing

.. toctree::
   :maxdepth: 2
   :caption: Modules

   modules/coordinate
   modules/encounter
   modules/satellite
   modules/schedule
   modules/time
   modules/units
