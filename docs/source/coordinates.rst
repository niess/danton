Coordinates systems
===================

Danton uses three system of coordinates, described below, in order to specify
the position and direction of Monte Carlo particles w.r.t. the
:doc:`Earth geometry <geometry>`.

.. topic:: Ellipsoid

   Danton coordinates systems all refer to a reference `ellipsoid`_. Depending
   on the selected Earth shape, the reference `ellipsoid`_ might or might not
   coincide with the `geoid`_ (see :numref:`tab-geoid`).


Geographic coordinates
----------------------

The Danton interface usually expresses positions using GPS-like `geodetic`_
coordinates (latitude, longitude, altitude). Directions are expressed using
`horizontal`_ coordinates (azimuth, elevation), clock-wise w.r.t. the geographic
north. These sets of coordinates are denoted as `Geographic` by Danton.

.. note::

   Altitudes are expressed w.r.t. the `ellipsoid`_, not w.r.t. the `geoid`_.
   Note that this is usually the case for GPS altitudes. As a result, a zero
   altitude value might differ from the sea-level.


Geocentric coordinates
----------------------

Internally, Danton uses Earth-Centered, Earth-Fixed (`ECEF`_) cartesian
coordinates (also refered to as `geocentric`). Conversions between geographic
and `ECEF`_ coordinates can be performed with the :py:func:`from_ecef
<danton.Geometry.from_ecef>` and :py:func:`to_ecef <danton.Geometry.to_ecef>`
methods of the :py:class:`Geometry <danton.Geometry>` class. Note that the
result depends on the `ellipsoid`_, thus on the selected `geoid`_.

.. note::

   Except for the :python:`"PREM"` ellipsoid (i.e. a spherical Earth), the
   geographic (radial) and geocentric (orthometric) altitudes differ.


Local coordinates
-----------------

In some cases, it is more convenient to use local cartesian coordinates, e.g.
when considering a box-bounded area of the Earth surface. In this case, Danton
uses Local-Tangent-Plane (`LTP`_) coordinates. Conversions between local and
geographic coordinates can be performed with the :py:func:`from_local
<danton.Box.from_local>` and :py:func:`to_local <danton.Box.to_local>` methods
of the :py:class:`Box <danton.Box>` class. Note that, as for geocentric
coordinates, the result depends on the `ellipsoid`_.

.. topic:: Orientation

   By default, Danton local coordinates are East-North-Upward (ENU) oriented,
   w.r.t. the geographic north. Optionally, a :py:attr:`declination
   <danton.Box.declination>` angle can be specified, clock-wise around the
   z-axis.


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _ECEF: https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system
.. _ellipsoid: https://en.wikipedia.org/wiki/Earth_ellipsoid
.. _geoid: https://en.wikipedia.org/wiki/Geoid
.. _geodetic: https://en.wikipedia.org/wiki/Geodetic_coordinates
.. _horizontal: https://en.wikipedia.org/wiki/Horizontal_coordinate_system
.. _LTP: https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates
