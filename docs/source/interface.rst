Python interface
================

.. autoclass:: danton.Box

   .. automethod:: from_local
   .. automethod:: inside
   .. automethod:: to_local

   .. autoattribute:: altitude
   .. autoattribute:: declination
   .. autoattribute:: geodesic
   .. autoattribute:: latitude
   .. autoattribute:: longitude
   .. autoattribute:: size
   .. autoattribute:: surface

----

.. autoclass:: danton.Geometry

   This class provides an interface to the Monte Carlo geometry, which may be
   configured through attributes. See the :doc:`geometry` section for more
   details.

   .. method:: __new__()

      Create a new Earth geometry.

      >>> geometry = danton.Geometry()

      By default, the Earth geometry is initialised according to the Preliminary
      Reference Earth Model (`PREM`_). See the attributes below for other
      options.

   .. autoattribute:: density

      Get or set the density of the topography layer, in kg/m\ :sup:`3`. For
      example, the following sets the ground density to 0.92 g/cm\ :sup:`3`.

      >>> geometry.density = 0.92E+03

   .. autoattribute:: geodesic

      Get or set the reference model for the sea level. By default, a spherical
      Earth is assumed, according to the `PREM`_. Possible values are listed
      in :numref:`tab-geodesic` below.

      .. _tab-geodesic:

      .. list-table:: Available geodesic models.
         :width: 75%
         :widths: auto
         :header-rows: 1

         * - Model
           - Description
         * - :python:`"EGM96"`
           - An undulated Earth, according to the `EGM96`_ model.
         * - :python:`"PREM"`
           - A spherical Earth, according to the `PREM`_ model.
         * - :python:`"WGS84"`
           - An elliptical Earth, according to the `WGS84`_ ellipsoid.

   .. autoattribute:: material

      Get or set the constitutive material of the topography layer, which is
      `standard rock`_ by default. For example, the following sets the ground
      composition to water (e.g. to represent an ice cover).

      >>> geometry.material = "Water"

   .. autoattribute:: ocean

      By default, the Earth is covered by a 3 km deep ocean, according to the
      `PREM`_. Setting this attribute to :python:`False` overrides the ocean by
      extending the upper crust layer up to the sea level.

   .. autoattribute:: topography

      Get or set the elevation model for describing the topography layer. By
      default, this is :python:`None`, i.e. there are no topography elements.
      A :py:class:`float` value specifies a uniform topography level, in m. For
      instance, the following sets the topography height to 1 km (all over the
      Earth).

      >>> geometry.topography = 1E+03

      Alternatively, a path to a folder containing tiles from a Global Digital
      Elevation Model (GDEM) can be provided (see e.g. `SRTMGL1`_). For
      instance, as

      >>> geometry.topography = "/path/to/gdem/"

----

.. autoclass:: danton.ParticlesGenerator

   .. automethod:: direction
   .. automethod:: energy
   .. automethod:: generate
   .. automethod:: inside
   .. automethod:: pid
   .. automethod:: position
   .. automethod:: powerlaw
   .. automethod:: solid_angle
   .. automethod:: target

----

.. autoclass:: danton.Physics

   .. autoattribute:: bremsstrahlung
   .. autoattribute:: dis
   .. autoattribute:: pair_production
   .. autoattribute:: pdf
   .. autoattribute:: photonuclear

----

.. autoclass:: danton.Random

   .. automethod:: uniform01
   .. autoattribute:: index
   .. autoattribute:: seed

----

.. autoclass:: danton.Simulation

   This class provides an interface to a Danton Monte Carlo simulation, which
   may be configured through attributes.

   .. method:: __new__()

      Create a new simulation interface, configured by default for the backward
      sampling of tau decays. For instance,

      >>> simulation = danton.Simulation()

   .. automethod:: box
   .. automethod:: particles

   .. automethod:: run

      The provided *particles* are transported through the Monte Carlo
      :py:attr:`geometry`. The returned object depends on the simulation
      :py:attr:`mode`, :py:attr:`record_steps` and :py:attr:`tau_decays`
      attributes. For example, in :python:`Backward` mode with tau decays (but
      no steps recording), a :external:py:class:`NamedTuple <typing.NamedTuple>`
      is returned, containing the sampled primaries, as well as the tau creation
      vertices, and the tau decay products (as
      :external:py:class:`numpy.ndarray`, each).

   .. autoattribute:: geometry

      Get, modify, or set the Monte Carlo :py:class:`Geometry`. For instance,
      the following sets an elliptical Earth, according to the `WGS84`_
      ellipsoid.

      >>> simulation.geometry.geodesic = "WGS84"

   .. autoattribute:: longitudinal

      By default, a full 3D Monte Carlo simulation is performed. When this
      attribute is set to :python:`True`, then deflections are neglected during
      collision (but not the resulting energy losses), leading to all particles
      propagating along the same line.

   .. autoattribute:: mode

      Must be one of :python:`"backward"` (default), :python:`"forward"` or
      :python:`"grammage"`.

   .. autoattribute:: physics
   .. autoattribute:: random

   .. autoattribute:: record_steps

      By default, only a few Monte Carlo states of interest are recorded (see
      e.g. the :py:class:`run` method). Setting this attribute to :python:`True`
      results in the full history of Monte Carlo steps being recorded.

   .. autoattribute:: tau_decays

      By default, Danton is configured to sample tau decays over a volume.
      Setting this flag to :python:`False` results in the sampling of flux
      events instead, through a boundary geodesic-like surface, defined by a
      constant altitude.

----

.. autofunction:: danton.particles

   This function returns a `structured numpy array <StructuredArray_>`_ with the
   given *shape*. Particles are initialised with default properties, if
   not overriden by specifying *kwargs*. For instance, the following creates an
   array of 100 particles (taus, by default) with an energy of
   1 EeV, starting from the Equator (by default), and going to the North.

   >>> particles = danton.particles(100, energy=1E+09)

   The data structure (:external:py:class:`numpy.dtype`) of a particle
   is the following (the corresponding physical units are also indicated).

   .. list-table:: Particles array structure.
      :width: 50%
      :widths: auto
      :header-rows: 1

      * - Field
        - Format
        - Units
      * - :python:`"pid"`
        - :python:`"i4"`
        - 
      * - :python:`"energy"`
        - :python:`"f8"`
        - GeV
      * - :python:`"latitude"`
        - :python:`"f8"`
        - deg
      * - :python:`"longitude"`
        - :python:`"f8"`
        - deg
      * - :python:`"altitude"`
        - :python:`"f8"`
        - m
      * - :python:`"azimuth"`
        - :python:`"f8"`
        - deg
      * - :python:`"elevation"`
        - :python:`"f8"`
        - deg

   .. topic:: Particle ID

      The type of a Monte Carlo particle (:python:`"pid"`) is encoded according
      to the Particle Data Group (PDG) `numbering scheme <PdgScheme_>`_.

.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _EGM96: https://cddis.nasa.gov/926/egm96/egm96.html
.. _PREM: https://en.wikipedia.org/wiki/Preliminary_reference_Earth_model
.. _PdgScheme: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
.. _SRTMGL1: https://lpdaac.usgs.gov/products/srtmgl1v003/
.. _standard rock: https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/standard_rock.html
.. _StructuredArray: https://numpy.org/doc/stable/user/basics.rec.html
.. _WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
