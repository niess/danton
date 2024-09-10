Geometry description
====================

Earth structure
---------------

Danton models the Earth's interior structure according to the Preliminary
Reference Earth Model (`PREM`_). That is, the Earth is composed of concentric
layers, ranging from the very core to the surface, with varying bulk densities
depending on the depth.

.. note::

   The legacy `PREM`_ Earth is covered by a 3 |nbsp| km deep ocean, without any
   land. This can be disabled however (see the :py:attr:`ocean
   <danton.Geometry.ocean>` attribute of the :py:class:`Geometry
   <danton.Geometry>` class).

.. topic:: Composition

   The `PREM`_ is derived from seismic data and does not specify the composition
   of the Earth. Thus, Danton assumes a uniform composition aligned with
   `standard rock`_, with the exception of the optional ocean layer.


Earth shape
-----------

Danton model's the sea level according to a `geoid`_. Three models are supported
(see the :py:attr:`geoid <danton.Geometry.geoid>` attribute of the
:py:class:`Geometry <danton.Geometry>` class).

The default is to consider a spherical `geoid`_ of radius 6371 |nbsp| km,
according to the `PREM`_. As a refinement, the `geoid`_ can instead be modelled
by an ellispoid, specifically the `WGS84`_ ellipsoid. An additional level of
accuracy is provided by the `EGM96`_ model (also included in Danton), whose
undulations are defined with reference to the `WGS84`_ ellipsoid.

.. topic:: Scaling

   The Earth interior structure is scaled radialy such that the ocean layer
   (i.e. :math:`x=1`, according to the `PREM`_) matches the `geoid`_. Note
   however that, except for a spherical geoid, this does not result in inner
   layers interfaces having a constant altitudes (or depth) w.r.t the reference
   ellipsoid.


Topography
----------

The Earth's surface can be modelled according to a Global Digital Elevation
Model (GDEM), such as `SRTMGL1`_. Alternatively, a constant topography elevation
can be specified for the entire Earth.

.. important::

   Danton always assumes that topography elevation values are given w.r.t. the
   `geoid`_.

.. topic:: Oceanic regions

   In the event that a topographic model is being specified, the `PREM`_ ocean
   layer is removed. However, negative topography elevations (e.g. from
   bathymetry data) indicate oceanic regions, ranging from the negative
   topography elevation values up to :math:`z=0`\, which corresponds to the
   `geoid`_.


Atmosphere
----------

Danton considers that the Earth's atmosphere adheres to the U.S. Standard
Atmosphere (`USSA`_) density profile (utilising Corsika's parameterization),
assuming a uniform composition. As for the Earth's interior, the atmosphere is
partitioned into layers, which are scaled to the `geoid`_.

.. topic:: Outer space

   The Monte Carlo geometry extends up to an altitude of 100 |nbsp| km
   (unscaled), which is considered by Danton to be the boundary with outer
   space.


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _EGM96: https://cddis.nasa.gov/926/egm96/egm96.html
.. _geoid: https://en.wikipedia.org/wiki/Geoid
.. _PREM: https://en.wikipedia.org/wiki/Preliminary_reference_Earth_model
.. _SRTMGL1: https://lpdaac.usgs.gov/products/srtmgl1v003/
.. _standard rock: https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/standard_rock.html
.. _USSA: _https://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
.. _WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
