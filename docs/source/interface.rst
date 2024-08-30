Python interface
================

.. autoclass:: danton.Box

   .. autoattribute:: altitude
   .. autoattribute:: declination
   .. autoattribute:: ellipsoid
   .. autoattribute:: latitude
   .. autoattribute:: longitude
   .. autoattribute:: size
   .. autoattribute:: surface

   .. automethod:: from_local
   .. automethod:: inside
   .. automethod:: to_local

----

.. autoclass:: danton.Geometry

   This class provides an interface to the Monte Carlo geometry, which may be
   configured through attributes. See the :doc:`geometry` section for more
   details.

   .. method:: __new__(**kwargs)

      Create a new Earth geometry.

      >>> geometry = danton.Geometry(geoid="EGM96")

      By default, the Earth geometry is initialised according to the Preliminary
      Reference Earth Model (`PREM`_). See the attributes below for other
      options.

   .. autoattribute:: density

      Get or set the density of the topography layer, in kg/m\ :sup:`3`. For
      example, the following sets the ground density to 0.92 g/cm\ :sup:`3`.

      >>> geometry.density = 0.92E+03

   .. autoattribute:: ellipsoid

      .. note::

         This attribute is read-only. It is defined by the :py:attr:`geoid`
         model (see :numref:`tab-geoid` below).

   .. autoattribute:: geoid

      Get or set the reference geoid for the sea level. By default, a spherical
      Earth is assumed, according to the `PREM`_. Possible values are listed
      in :numref:`tab-geoid` below.

      .. _tab-geoid:

      .. list-table:: Available geoid models.
         :width: 75%
         :widths: auto
         :header-rows: 1

         * - Model
           - Description
           - Ellipsoid
         * - :python:`"EGM96"`
           - An undulated Earth, according to the `EGM96`_ model.
           - :python:`"WGS84"`
         * - :python:`"PREM"`
           - A spherical Earth, according to the `PREM`_ model.
           - :python:`"PREM"`
         * - :python:`"WGS84"`
           - An elliptical Earth, according to the `WGS84`_ ellipsoid.
           - :python:`"WGS84"`

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

   .. automethod:: from_ecef

      The *position* and *direction* arguments must be arrays-like containing
      cartesian coordinates in Earth-Centered Earth-Fixed (ECEF) frame. For
      instance,

      >>> geodetic_position = geometry.from_ecef([6378137, 0, 0])

   .. automethod:: geoid_undulation

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a structured
      :external:py:class:`numpy.ndarray` containing geodetic coordinates
      (:python:`"latitude"`, :python:`"longitude"`).

      Alternativelly, geodetic coordinates can be specified as *kwargs*. For
      instance,

      >>> undulation = geometry.geoid_undulation(latitude=45, longitude=3)

   .. automethod:: locate

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      If *resolve* is :python:`False` (the default), then atmospheric layers are
      undifferentiated.

   .. automethod:: to_ecef

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a structured
      :external:py:class:`numpy.ndarray` containing geodetic coordinates
      (:python:`"latitude"`, :python:`"longitude"`, :python:`"altitude"`), and
      optionally horizontal ones (:python:`"azimuth"`, :python:`"elevation"`).

      Alternativelly, geodetic (and horizontal) coordinates can be specified as
      *kwargs*. For instance,

      >>> ecef_position = geometry.to_ecef(latitude=45, longitude=3, altitude=0)

   .. automethod:: topography_elevation

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a structured
      :external:py:class:`numpy.ndarray` containing geodetic coordinates
      (:python:`"latitude"`, :python:`"longitude"`).

      Alternativelly, geodetic coordinates can be specified as *kwargs*. For
      instance,

      >>> z = geometry.topography_elevation(latitude=45, longitude=3)

      The optional *reference* argument specifies the reference surface for
      elevation values. Possible values are :python:`"ellipsoid"` (default) or
      :python:`"geoid"`.

   .. automethod:: trace

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a structured
      :external:py:class:`numpy.ndarray` containing start position
      (:python:`"latitude"`, :python:`"longitude"`, :python:`"altitude"`) and
      tracing direction (:python:`"azimuth"`, :python:`"elevation"`).

      Alternativelly, position and direction can be specified as *kwargs*. For
      instance,

      >>> trace = geometry.trace(
      ...     latitude = 45, longitude = 3, altitude = 0,
      ...     azimuth = 0, elevation = 90
      ... )

      The optional *backward* argument reverses the tracing direction, if set to
      :python:`True`.

      The optional *limit* argument specifies a maximum distance, in m, for the
      tracing.

      If *resolve* is :python:`False` (the default), then atmospheric layers are
      undifferentiated.


   .. automethod:: translate

      Translate geodetic coordinates by *distance* (in m) along a straight line.
      If a negative distance is provided, then the translation direction is
      reversed.

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The positional *array* argument, if specified, must be a structured
      :external:py:class:`numpy.ndarray` containing the initial position
      (:python:`"latitude"`, :python:`"longitude"`, :python:`"altitude"`) and
      the translation direction (:python:`"azimuth"`, :python:`"elevation"`).

      Alternativelly, the initial position and the translation direction can be
      specified as *kwargs*. For instance,

      >>> coordinates = geometry.translate(
      ...     1000,
      ...     latitude = 45, longitude = 3, altitude = 0,
      ...     azimuth = 0, elevation = 90
      ... )

      When an *array* argument is provided, set *copy* to :python:`False` in
      order to translate coordinates in-place.

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

   .. autoattribute:: index
   .. autoattribute:: seed

   .. automethod:: uniform01

----

.. autoclass:: danton.Simulation

   This class provides an interface to a Danton Monte Carlo simulation, which
   may be configured through attributes.

   .. method:: __new__()

      Create a new simulation interface, configured by default for the backward
      sampling of tau decays. For instance,

      >>> simulation = danton.Simulation()

   .. autoattribute:: geometry

      Get, modify, or set the Monte Carlo :py:class:`Geometry`. For instance,
      the following sets an elliptical Earth, according to the `WGS84`_
      ellipsoid.

      >>> simulation.geometry.geoid = "WGS84"

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
      events instead, through a boundary surface defined by a constant altitude
      (w.r.t. the reference ellipsoid).

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
