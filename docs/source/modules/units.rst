Units
=====

.. contents:: Contents
   :depth: 1
   :local:

Units Overview
--------------

Celest provides the ability to handle different units of measurement and is accomplished through the use of units an
quantities. The following sections describe the different units available, how to use them, and how to convert between
them via the :class:`Quantity` class.

Handling Units
--------------

Initialized units are accessed from the units module which can be imported by the following:

.. code-block:: python

   from celest import units as u

Primary Units
^^^^^^^^^^^^^

Primary units can then be accessed using a dot notation.

.. code-block:: python

   u.m  # Meter measure.
   u.s  # Second measure.
   u.deg # Degree measure.

The following table lists all primary units, their code call sign, and their measure type.

.. list-table:: Available Units
   :widths: 30, 30, 30
   :header-rows: 1

   * - Unit
     - Code Label
     - Type
   * - Millimeter
     - `u.mm`
     - Length
   * - Centimeter
     - `u.cm`
     - Length
   * - Meter
     - `u.m`
     - Length
   * - Kilometer
     - `u.km`
     - Length
   * - Inch
     - `u.inch`
     - Length
   * - Feet
     - `u.ft`
     - Length
   * - Yard
     - `u.yd`
     - Length
   * - Mile
     - `u.mi`
     - Length
   * - Second
     - `u.s`
     - Time
   * - Minute
     - `u.min`
     - Time
   * - Hour
     - `u.hr`
     - Time
   * - Day
     - `u.dy`
     - Time
   * - JD2000
     - `u.jd2000`
     - Date
   * - Degree
     - `u.deg`
     - Angle
   * - Radian
     - `u.rad`
     - Angle
   * - Seconds of Arc
     - `u.arcsec`
     - Angle
   * - Minutes of Arc
     - `u.arcmin`
     - Angle
   * - Hour Angle
     - `u.hourangle`
     - Angle

Primary units are instances of the :class:`Unit` class and can be used to create new units.

.. autoclass:: celest.units.core.Unit
   :members:
   :inherited-members:

Compound Units
^^^^^^^^^^^^^^

Primary units can be combined with others through multiplication or division to create compound units.

.. code-block:: python

   velocity_unit = u.m / u.s  # Meter per second.
   acceleration_unit = u.m / u.s ^ 2  # Meter per second squared.

Compound units are instances of the :class:`CompoundUnit` class and can be used to create new units.

.. autoclass:: celest.units.core.CompoundUnit
   :members:
   :inherited-members:

Passing Units into Celest
^^^^^^^^^^^^^^^^^^^^^^^^^

Units are typically passed into Celest's class constructors, methods, and functions as a parameter to specify the data
measure. Internally, the data and units get stored together in a :class:`Quantity` object that may be returned from a
class property, method, or function. The :class:`Quantity` object deals with the conversions between units and is
discussed the :ref:`Measures with Units` section.

Measures with Units
-------------------

To allow for unit conversions, Celest uses the :class:`Quantity` class to store data and units together. This class is
typically initialized internally by a class constructor, method, or function and returned to users. The :class:`Quantity`
class documentation follows.

.. TODO: Detail what constraints apply to the data such that unit conversions can occur (e.g. an addition dunder).

.. autoclass:: celest.units.quantity.Quantity
   :members:
   :noindex:
