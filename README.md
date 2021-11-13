# Celest 0.2.0

Celest is a satellite dynamics and mission planning library designed for the University of Toronto Aerospace Team which is applying the program to a hyperspectral imaging CubeSat mission. The library aims to provide the necessary tools to plan satellite-ground interactions such as ground target tracking, and calculating imaging and transmission opportunities.

* **Source Code:** https://github.com/JaiWillems/Celest
* **Bug Report or Feature Request:** https://github.com/JaiWillems/Celest/issues

## Installation
Celest can be installed from PyPI using the following command:

```terminal
pip install celest
```

## What can the program do?

The program only requires an input data set for a satellite that contains the positions in either a gcrs, itrs, or geographical frame and the assosciated Julian times. The program can then be used to calculate various position and time representations, satellite to ground encounters, and transformations necessary for ground location tracking.

## What can be expected of the project?

In future releases, a high order c-implemented orbit propagator will be incorporated into the library to allow for a complete mission planning workflow from inputting orbital parameters to calculating satellite-to-ground encounter opportunities.

We also intend to incorporate a window scheduling system that will be used to develop an optimal satellite itinerary for satellite mission planning.

Other features that are being considered include a terminal interface, higher fidelity conversion calculations, window/data visuals, data loaders, and more robust ground tracking by incorporating forward motion compensation.
