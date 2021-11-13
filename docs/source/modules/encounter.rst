Encounter Module
================

GroundPosition Class
--------------------

The `GroundPosition` class is used to represent a ground location that can be
passed into window generation functionality to produce viable ground-satellite
encounter opportunities. Currently, the `GroundPosition` class is fairly
preliminary but will be expanded in the future to accommodate more specific and
complicated ellipsoid models based on the input geographical coordinates.

.. autoclass:: celest.encounter.GroundPosition
   :members:
   :show-inheritance:
   :noindex:

Window Generation
-----------------

The window generation function takes produces all viable encounter
opportunities of a specified type for a given satellite and ground location.

.. autofunction:: celest.encounter.windows.generate
   :noindex:

Windows Class
-------------

The `Windows` class is the returned data structure from the `generate` function
which holds all encounter opportunities and provides some basic functionality
to gain insight or access to the data.

.. autoclass:: celest.encounter._window_handling.Windows
   :members:
   :show-inheritance:
   :noindex:
