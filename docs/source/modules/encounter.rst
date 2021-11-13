Encounter Module
================

GroundPosition Class
--------------------

The :class:`GroundPosition` class is used to represent a ground location that
can be passed into window generation functionality to produce viable
ground-satellite encounter opportunities. Currently, the
:class:`GroundPosition` class is fairly preliminary but will be expanded in the
future to accommodate more specific and complicated ellipsoid models based on
the input geographical coordinates.

.. autoclass:: celest.encounter.GroundPosition
   :members:
   :show-inheritance:
   :noindex:

Window Generation
-----------------

The window generation function produces all viable encounter opportunities of a
specified type for a given satellite and ground location.

.. autofunction:: celest.encounter.windows.generate
   :noindex:

Window Handling
---------------

The :class:`Windows` class is the returned data structure from the
:py:func:`windows.generate` function which holds all encounter opportunities as
:class:`Window` objects and provides basic functionality to gain insight and
access to the data. The following details these two classes.

Windows
~~~~~~~

The :class:`Windows` class is interfaced through a series of interactive
methods. The window data can be interfaced in two ways. The first is through
the :class:`Windows.windows` attribute which holds the :class:`Window`
objects in a Pandas :py:func:`Series` data structure indexed by the window
start time. Window data can also be accessed by iterating over the class.

.. code-block::

   # Window accessed by indexing by the start time 2400000.5
   window_object = windows_object.windows[2400000.5]

   # Window accessed by iterating over the class.
   for window_object in windows_object:

      # window_object window data can be accessed using attributes.
      lat, lon = window_object.coor
      start = window_object.start
      end = window_object.start
   

.. autoclass:: celest.encounter._window_handling.Windows
   :members:
   :show-inheritance:
   :noindex:

Window
~~~~~~

The :class:`Window` class is the data structure that holds all the information
regarding a specific encounter opportunity. The important information is held
in a series of attributes.

.. autoclass:: celest.encounter._window_handling.Window
   :members:
   :show-inheritance:
   :noindex:
