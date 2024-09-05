Geometry description
====================

Earth structure
---------------

Danton Model's the Earth interior structure according to the Preliminary
Reference Earth Model (`PREM`_). That is, the Earth is composed of concentric
layers, ranging from the very core to the surface, with varying bulk densities
depending on the depth.

.. note::

   The legacy `PREM`_ Earth is covered by a 3 |nbsp| km deep ocean, without any
   land. This can be disabled however (see the :py:attr:`ocean
   <danton.Geometry.ocean>` attribute of the :py:class:`Geometry
   <danton.Geometry>` class).

.. topic:: Composition

   The PREM originates from seismic data. It does not specify the Earth's
   composition, but only seismic related profiles (the like, bulk density).
   Thus, Danton currently assumes a uniform composition corresponding to
   `standard rock`_, except for the optional ocean layer, for which pure water
   is assumed (but with PREM's / salted water density).


Earth shape
-----------

Danton model's the sea level according to a `geoid`_. Three models are supported
(see the :py:attr:`geoid <danton.Geometry.geoid>` attribute of the
:py:class:`Geometry <danton.Geometry>` class).

The default is to consider a spherical `geoid`_ of radius 6371 |nbsp| km,
according to the `PREM`_. As a refinement, the `geoid`_ can instead be modelled
by an ellispoid, the `WGS84`_ one namely. An additional level of accuracy is
provided by the `EGM96`_ model (also included in Danton), whose undulations are
defined w.r.t. the `WGS84`_ ellipsoid.

.. topic:: Scaling

   The Earth interior structure is scaled radialy such that the ocean layer
   (i.e. :math:`x=1`, according to the `PREM`_) matches the `geoid`_. Note
   however that, except for a spherical geoid, this does not result in inner
   layers interfaces having a constant altitudes (or depth) w.r.t the reference
   ellipsoid.


Topography
----------

Optionally, the Earth surface can be modelled by a Global Digital Elevation
Model (GDEM), the like e.g. `SRTMGL1`_. Alternativelly, a constant topography
elevation, all over the Earth, can be specified as well.

.. important::

   Danton always assumes that topography elevation values are given w.r.t. the
   `geoid`_.

.. topic:: Oceanic regions

   When specifying a topographic model, the `PREM`_ ocean layer is removed.
   However, negative topography elevations (e.g. from bathymetry data) indicate
   oceanic regions, ranging from the negative topography elevation values up to
   :math:`z=0`\, which corresponds to the geoid.


Atmosphere
----------

Danton considers that the Earth atmopshere follows the U.S. Standard (`USSA`_)
density profile (using Corsika's parameterisation), with a uniform composition.
As for the Earth interior, the atmosphere is divided in layers, which are scaled
according to the `geoid`_.

.. topic:: Outer space

   The Monte Carlo geometry extends up to an altitude of 100 |nbsp| km
   (unscaled), which is considered by Danton to be the boundary with the outer
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
