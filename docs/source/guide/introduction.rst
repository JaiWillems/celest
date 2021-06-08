Introduction
============

Celest is a Python library designed for satellite position representations, coordinate conversions, and encounter planning.

What seperates Celest from other astronomical libraries are its computational efficiencies acheived by model simplifications, leveraging Numpy's small overhead, and vectorizing inputs where ever possible. With these considerations in mind, the Celest library does not compromize accuracy. For example, Celest's ECI to horizontal conversion have 100% similarity with Astropy with only a 0.0175% tolerance.

Motivation
**********

Originally designed for the University of Toronto Aerospace Team's Space Systems division, Celest was a neccessary addition to the teams kit out of a desire for code flexibility but most fundamentally, efficiency.

Still vastly customized for a cube-satellite remote sensing mission, the horizons of Celest are focussing on broadening the scope of the software to increase applications and reduce the need to use various libraries in parallel.

What can Celest Do?
*******************

#. Data Abstraction
    * Celest provides various objects to consolidate information and simplify conversions.
#. Conversions
    * Celest can convert between ECI, ECEF, and horizontal coordinate frames,
#. Encounter Planning
    * Apply sun and angular constraints to generate encounter windows,
    * Employ an adaptive interpolation method to increase window precision and reduce runtime.
