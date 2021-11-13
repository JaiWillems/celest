Quickstart
==========

Installation using pip
----------------------

If you use `pip`, you can install Celest with:

.. code-block::

   pip install celest

Getting familiar with the library
---------------------------------

Celest is broken primarily into two modules: (1) :mod:`celest.satellite` module
contains all functions related to satellite dynamics such as
coordinate and time conversions. (2) :mod:`celest.encounter` module contains
all functions related to the encounters between a satellite and ground location.

A typical workflow for using Celest is to begin by importing satellite time and
position data as NumPy arrays which can then be formed into :class:`Time` and
:class:`Coordinate` objects which are then used to initialize a
:class:`Satellite` object. The :class:`Satellite` object can then be coupled
with encounter information such as a ground location and encounter parameters
to determine encounter opportunities.

.. code-block::

   import numpy as np
   from celest.satellite import Time, Coordinate, Satellite
   from celest.encounter import GroundPosition, windows

   # Load the data.
   julian = np.loadtxt('julian.txt')
   gcrs = np.loadtxt('gcrs.txt')

   # Initialize satellite representation.
   time = Time(julian, offset=0)
   position = Coordinate(gcrs, frame='gcrs', time=time)
   satellite = Satellite(position=position)

   # Define ground position.
   toronto = GroundPosition(latitude=43.65, longitude=-79.38)

   # Generate ground location windows.
   toronto_IMG_windows = windows.generate(satellite=satellite, location=toronto, enc="image", ang=30)
   toronto_GL_windows = windows.generate(satellite=satellite, location=toronto, enc="data link", ang=10)

   # Save satellite encounter windows.
   toronto_IMG_windows.save(fname="toronto_IMG_windows.csv", delimiter=",")
   toronto_DL_windows.save(fname="toronto_DL_windows.csv", delimiter=",")
