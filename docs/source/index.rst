.. Celest documentation master file, created by
   sphinx-quickstart on Mon Jun  7 18:56:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Celest Documentation
====================

Celest is a satellite mission planning program designed for the University of
Toronto Aerospace Team which is applying the program to a hyperspectral
imaging CubeSat mission. The library aims to provide the necessary tools to
plan satellite-ground interactions such as ground target tracking, and
calculating imaging and ground transmission opportunities.

Being subjected to a design team workflow, Celest is consistently under
development for improvements in efficiency, functionality, and fidelity. It is 
then expected that the library will be updated regularly.

What can the program do?
------------------------

The program only requires an input data set for a satellite that contains the
positions in either a gcrs, itrs, or geographical coordinates and the
assosciated Julian times. The program can then be used to calculate various
position and time representations, satellite to ground encounters, and
transformations necessary for ground location tracking.


What can be expected of the project?
------------------------------------

In future releases, a high order c-implemented orbit propagator will be
incorporated into the library to allow for a complete mission planning workflow
from inputting orbital parameters to calculating satellite-to-ground encounter
opportunities.

We also intend to incorporate a window scheduling system that will be used to
develop an optimal satellite itinerary for satellite mission planning.

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
