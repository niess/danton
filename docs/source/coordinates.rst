Coordinates systems
===================

Danton uses three distinct systems of coordinates, which are outlined below, to
define the position and direction of Monte Carlo particles in relation to the
:doc:`Earth's geometry <geometry>`.

.. topic:: Ellipsoid

   Danton coordinate systems are based on a reference `ellipsoid`_. Depending on
   the selected Earth shape, the reference `ellipsoid`_ may or may not align
   with the `geoid`_ (see see :numref:`tab-geoid`).


Geographic coordinates
----------------------

The Danton interface primarily uses GPS-like `geodetic`_ coordinates (latitude,
longitude, altitude) to express positions. Directions are expressed using
`horizontal`_ coordinates (azimuth, elevation), clockwise w.r.t. the geographic
north. These sets of coordinates are collectively referred to as `Geographic` by
Danton.

.. note::

   Altitudes are expressed w.r.t. the `ellipsoid`_ rather than the `geoid`_.
   This is a standard practice for GPS altitudes. Consequently, a zero altitude
   value may differ from the sea-level.


Geocentric coordinates
----------------------

Internally, Danton uses Earth-Centred, Earth-Fixed (`ECEF`_) Cartesian
coordinates, also known as geocentric. The  :py:func:`from_ecef
<danton.Geometry.from_ecef>` and :py:func:`to_ecef <danton.Geometry.to_ecef>`
methods of the :py:class:`Geometry <danton.Geometry>` class facilitate
conversions between geographic and `ECEF`_ coordinates. It should be noted that
the result is dependent on the `ellipsoid`_, which in turn depends on the
selected `geoid`_.

.. note::

   With the exception of the :python:`"PREM81"` ellipsoid (i.e. a spherical
   Earth), there is a discrepancy between the geographic (orthometric) and
   geocentric (radial) altitudes.


Local coordinates
-----------------

In some cases, it is more convenient to use local cartesian coordinates, e.g.
when considering a box-bounded area of the Earth's surface. In this case, Danton
uses Local-Tangent-Plane (`LTP`_) coordinates. Conversions between local and
geographic coordinates can be performed with the :py:func:`from_local
<danton.Box.from_local>` and :py:func:`to_local <danton.Box.to_local>` methods
of the :py:class:`Box <danton.Box>` class. Note that, as with geocentric
coordinates, the result depends on the `ellipsoid`_.

.. topic:: Orientation

   By default, Danton local coordinates are East-North-Upward (ENU) oriented,
   w.r.t. the geographic north. Optionally, a :py:attr:`declination
   <danton.Box.declination>` angle can be specified, in a clock-wise direction
   around the z-axis.


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
