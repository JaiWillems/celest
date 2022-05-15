Encounter Module
================

GroundPosition Class
--------------------

The :class:`GroundPosition` class is used to represent a ground location that
can be passed into window generation functionality to produce viable
ground-satellite encounter opportunities. Currently, the :class:`GroundPosition`
class is somewhat preliminary but will be expanded in the future to accommodate
more complicated ellipsoid models based on the input geographical coordinates.

.. autoclass:: celest.encounter.GroundPosition
   :members:
   :show-inheritance:
   :noindex:

Window Generation
-----------------

The window generation function determines all visible time windows for a
satellite and ground-location pair where the satellites elevation is
greater than the visibility threshold.

.. autofunction:: celest.encounter.windows.generate_vtw
   :noindex:

Window Handling
---------------

The :class:`VTWHandler` class is the data structure returned from the
:py:func:`windows.generate_vtw` function which holds all visible time windows
as :class:`VTW` objects and provides an interface to access the data.
Similarly, the :class:`OWHandler` class is the data structure returned from
the scheduling workflow which hols all observation windows as :class:`OW`
objects. The following detail these various classes.

VTWHandler
~~~~~~~~~~

The window data held within the :class:`VTWHandler` class can be interfaced in
three ways: indexing, using access methods, and iterating.

To provide meaningful interactions with the window data, the
:class:`VTWHandler` is indexed by Julian start times in the J2000 epoch.
Although, the start times are typically unknown to the user. As a result, a
unique indexing scheme was adopted by the class to provide data access that
might be convenient in context. The window associated with the closest start
time is returned when indexing using an integer or float value. Similarly,
when indexing with a tuple of values, the unique windows related to the closest
start time are returned. When two indices are both closest to the start time of
a window, the window is only returned once. All windows with start times
falling within the slice (inclusive) are returned if the index is a slice.

This indexing scheme allows users to easily access the data in the
:class:`VTWHandler` class without prior knowledge of start times. It also
allows the user to determine the encounter closest to the desired time for
scheduling purposes.

The window data can also be interfaced with using the :py:func:`get_window` and
:py:func:`get_windows_in_range` methods. The first returns the window with the
closes start time. The latter returns all windows with start times falling
within the input range.

Lastly, window data can be accessed by iterating over the class.

Examples of object interactions are seen in the following example.

.. code-block::

   # Assuming we have a Windows object with three windows with start times
   # 2405795.5, 2405796.5, and 2405797.5.

   # Window accessed by indexing by the start time 2405796.5.
   window_object = windows_object[2405796.5]
   # window_object now contains the window starting at 2405796.5.

   # Window accessed by indexing with a tuple.
   window_object = windows_object[2405796.5, 2405796.7]
   # window_object is a list containing windows with the 2405795.5, 2405796.5 start times.

   # Window accessed by indexing with a slice.
   window_object = windows_object[2405796.5:2405796.7]
   # window_object is a list containing windows with the 2405795.5, 2405796.5 start times.

   # Window accessed through the get_window method by the start time 2405796.5.
   window_object = windows_object.get_window(2405796.5)
   # window_object now contains the window starting at 2405796.5.

   # Window accessed through the get_windows_in_range method using a start and end time.
   window_object = windows_object.get_windows_in_range(2405796.5, 2405796.7)
   # window_object is a list containing windows with the 2405795.5, 2405796.5 start times.

   # Window accessed by iterating over the class.
   for window_object in windows_object:

      # window_object window data can be accessed using attributes.
      lat, lon = window_object.coor
      start = window_object.start
      end = window_object.start
   

.. autoclass:: celest.encounter._window_handling.VTWHandler
   :members:
   :show-inheritance:
   :noindex:

VTW
~~~

The :class:`VTW` class is the data structure that holds all the information
regarding a specific encounter opportunity. The important information is held
in a series of attributes.

.. autoclass:: celest.encounter._window_handling.VTW
   :members:
   :show-inheritance:
   :noindex:

OWHandler
~~~~~~~~~

The :class:`OWHandler` is retuned from the :class:`Schedule` when the
:py:func:`generate` method is called. It holds the :class:`OW` objects which
represent the scheduled requests. The :class:`OWHandler` class can be
interfaced in a similar way to the :class:`VTWHandler` class.

.. autoclass:: celest.encounter._window_handling.OWHandler
   :members:
   :show-inheritance:
   :noindex:

OW
~~

The :class:`OW` class is the data structure that holds all the information
regarding a scheduled request. The important information is held
in a series of attributes.

.. autoclass:: celest.encounter._window_handling.OW
   :members:
   :show-inheritance:
   :noindex:
