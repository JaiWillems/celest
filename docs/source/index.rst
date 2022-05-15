.. Celest documentation master file, created by
   sphinx-quickstart on Mon Jun  7 18:56:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Celest Documentation
====================

Celest is a satellite mission planning software designed for remote sensing
CubeSat mission proifles. The library aims to provide the necessary tools for
celestial orbital conversions and satellite-to-ground interactions such as
determining and scheduling imaging/transmission opportunities.

What can the program do?
------------------------

Currently, the program requires satellite position and velocity data in either
the gcrs or itrs frames and the corresponding times. Celest can then calculate
different position representations, time representations, and can determine and
schedule satellite-to-ground encounters.

What can be expected of the project?
------------------------------------

In future releases, a high order c-implemented orbit propagator will be
incorporated into the library to allow for a complete mission planning workflow
from inputting orbital parameters to calculating satellite-to-ground encounter
opportunities.

Other features that are being considered include a terminal interface, higher
fidelity conversion calculations, window/data visuals, data loaders, and more
robust ground tracking by incorporating forward motion compensation.

.. toctree::
   :maxdepth: 2
   :caption: Guide

   guide/quickstart
   guide/tutorials
   guide/contributing

.. toctree::
   :maxdepth: 2
   :caption: Modules

   modules/satellite
   modules/encounter
   modules/schedule
