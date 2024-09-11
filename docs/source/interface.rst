Python interface
================

.. autoclass:: danton.Box

   This class defines a bounding box for an area of interest in the Monte Carlo
   geometry. It also allows for transformation between :ref:`geographic
   <coordinates:Geographic coordinates>` and :ref:`local <coordinates:Local
   coordinates>` coordinates.

   .. method:: __new__(size=None, **kwargs)

      Create a new bounding-box.

      The *size* argument defines the dimensions of the box along the local x, y
      and z axes. It can be a :py:class:`float` (cube), a size two array (square
      base) or a size three array. For examples, the following creates a 1
      |nbsp| km high box with a 10x10 |nbsp| km\ :sup:`2` base.

      >>> box = danton.Box([1E+04, 1E+03], latitude=45, longitude=3)

      Refer to the attributes below for the possible values of the optional
      *kwargs*.

   .. autoattribute:: altitude

      The altitude is expressed in m, w.r.t. the `ellipsoid`_ (see the
      :doc:`coordinates` section).

   .. autoattribute:: declination

      The declination angle is expressed in deg, w.r.t. the geographic north
      (see the :ref:`coordinates:Local coordinates` section).

   .. autoattribute:: ellipsoid

      See :numref:`tab-geoid` below for possible values.

   .. autoattribute:: latitude
   .. autoattribute:: longitude
   .. autoattribute:: size
   .. autoattribute:: surface_area
   .. autoattribute:: volume

   .. automethod:: from_local

      The *position* and *direction* arguments must be arrays-like containing
      :ref:`local <coordinates:Local coordinates>` Cartesian coordinates in the
      box frame. For instance,

      >>> coordinates = box.from_local([0, 0, 100])

   .. automethod:: inside

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      This method returns a boolean (or a :external:py:class:`ndarray
      <numpy.ndarray>` of booleans), where :python:`True` value(s) indicate that
      coordinates are inside the box.

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates (:python:`"latitude"`, :python:`"longitude"`,
      :python:`"altitude"`).

      Alternativelly, :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates can be specified as *kwargs*. For instance,

      >>> inside = box.inside(latitude=45, longitude=3, altitude=0)

   .. automethod:: projected_area

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates (:python:`"azimuth"`, :python:`"elevation"`).

      Alternativelly, :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates can be specified as *kwargs*. For instance,

      >>> area = box.projected_area(azimuth=30, elevation=45)

   .. automethod:: to_local

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates (:python:`"latitude"`, :python:`"longitude"`,
      :python:`"altitude"`, and optionally :python:`"azimuth"`,
      :python:`"elevation"`).

      Alternativelly, :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates can be specified as *kwargs*. For instance,

      >>> position = box.to_local(latitude=45, longitude=3, altitude=0)

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

   .. autoattribute:: ocean

      By default, the Earth is covered by a 3 km deep ocean, according to the
      `PREM`_. Setting this attribute to :python:`False` overrides the ocean by
      extending the upper crust layer.

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

   .. autoattribute:: topography_density

      Get or set the density of the topography layer, in kg/m\ :sup:`3`. For
      example, the following sets the ground density to 0.92 g/cm\ :sup:`3`.

      >>> geometry.density = 0.92E+03 # kg/m3

   .. autoattribute:: topography_material

      Get or set the constitutive material of the topography layer, which is
      :python:`"Rock"` by default. For instance, the following sets the ground
      composition to water (to represent an ice cover, for example).

      >>> geometry.material = "Water"

   .. automethod:: from_ecef

      The *position* and *direction* arguments must be arrays-like containing
      :ref:`geocentric <coordinates:Geocentric coordinates>` Cartesian
      coordinates in Earth-Centered Earth-Fixed (`ECEF`_) frame. For instance,

      >>> coordinates = geometry.from_ecef([6378137, 0, 0])

   .. automethod:: geoid_undulation

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing :ref:`geographic <coordinates:Geographic coordinates>`
      (:python:`"latitude"`, :python:`"longitude"`).

      Alternativelly, :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates can be specified as *kwargs*. For instance,

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

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates (:python:`"latitude"`, :python:`"longitude"`,
      :python:`"altitude"`, and optionally :python:`"azimuth"`,
      :python:`"elevation"`).

      Alternativelly, :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates can be specified as *kwargs*. For instance,

      >>> position = geometry.to_ecef(latitude=45, longitude=3, altitude=0)

   .. automethod:: topography_elevation

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates (:python:`"latitude"`, :python:`"longitude"`).

      Alternativelly, :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates can be specified as *kwargs*. For instance,

      >>> z = geometry.topography_elevation(latitude=45, longitude=3)

      The optional *reference* argument specifies the reference surface for
      elevation values. Possible values are :python:`"ellipsoid"` (default) or
      :python:`"geoid"`.

   .. automethod:: trace

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing start position (:python:`"latitude"`, :python:`"longitude"`,
      :python:`"altitude"`) and tracing direction (:python:`"azimuth"`,
      :python:`"elevation"`).

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

      Translate :ref:`geographic <coordinates:Geographic coordinates>`
      coordinates by *distance* (in m) along a straight line. If a negative
      distance is provided, then the translation direction is reversed.

      .. note::

         The positional *array* argument and keyword only (*kwargs*) arguments
         are mutually exclusive.

      The positional *array* argument, if specified, must be a `structured
      <StructuredArray_>`_ :external:py:class:`ndarray <numpy.ndarray>`
      containing the initial position (:python:`"latitude"`,
      :python:`"longitude"`, :python:`"altitude"`) and the translation direction
      (:python:`"azimuth"`, :python:`"elevation"`).

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

.. autoclass:: danton.Materials

   This class provides an interface to a set of Monte Carlo materials, which are
   specified by a Materials Description File (MDF), in `TOML`_ format. Refer to
   the :doc:`materials` section for further details.

   .. method:: __new__(path=None, \)

      Load a set of Monte Carlo materials from a MDF. For instance,

      >>> materials = danton.Materials("examples/rocks.toml")

      If *path* is :python:`None`, Danton's default materials are loaded.

   .. method:: __getitem__(name, \)

      Return the properties of a material given its *name*.

----

.. autoclass:: danton.ParticlesGenerator

   This class provides a utility for the generation of Monte Carlo particles
   from configurable distributions. This tool is typically used to seed the
   Monte Carlo simulation with an initial set of particles.

   Once initialised, the generator can be further configured using the methods
   detailed below (i.e. following a `builder`_ pattern). The
   :py:func:`generate` method then triggers the actual sampling of Monte Carlo
   particles. As an example, the following generates N particles uniformly
   within a solid angle (spanning elevation values between -5 and +5 |nbsp|
   deg), and with a :math:`1/E^2` power-law energy distribution (between 1
   |nbsp| PeV and 1 |nbsp| EeV).

   >>> particles = generator                    \
   ...     .solid_angle(elevation=[-5, 5])      \
   ...     .powerlaw(1E+06, 1E+12, exponent=-2) \
   ...     .generate(N)

   .. method:: __new__(*, geometry=None, random=None, weight=True)

      Create a new particles generator.

      The *weight* argument specifies whether the particles should be weighted
      by the inverse of their generation likelihoods (:math:`\omega = 1 /
      \text{pdf}(\text{S})`, for a Monte Carlo state :math:`\text{S}`) or not.
      Note that this can be overridden by individual distributions (using the
      :python:`weight` flag of other methods).

   .. automethod:: direction
   .. automethod:: energy

   .. automethod:: generate

      The *shape* argument defines the number of particles requested (as a
      :external:py:class:`ndarray <numpy.ndarray>` shape).

      The outcome of this method is dependent on the generator configuration. If
      rejection sampling techniques are employed, then the actual sample size is
      returned along with the selected particles. Otherwise, only the generated
      particles are returned.

   .. automethod:: inside

      .. note:: This method utilises rejection sampling techniques.

      This method is best used in conjunction with a backward simulation, as it
      generates potential tau decay vertices within the specified *box* volume,
      contingent on being located inside the :ref:`atmosphere
      <geometry:Atmosphere>`.

      By default, the candidate vertices are selected based on their likelihood
      of originating from the ground (including any topography), according to
      the tau lifetime. It is possible to disable this selection by setting the
      *limit* argument to :python:`False`.

      When setting the *limit* argument to a :py:class:`float`, the provided
      value controls the maximum distance (in multiples of the tau decay length)
      of candidate vertices to the ground. The default limit is 3 times the tau
      decay length (i.e., :python:`limit=3`).

   .. automethod:: pid

      Monte Carlo particles are indentified by their Particle ID (PID), which
      follows the Particle Data Group (`PDG <PdgScheme_>`_) numbering scheme.

   .. automethod:: position

   .. automethod:: powerlaw

      The *energy_min* and *energy_max* arguments define the support of the
      power-law, as an interval.

      The default setting is a :math:`1 / E` power-lawx, corresponding to
      :python:`exponent=-1`. Note that setting the exponent value to zero
      results in a uniform distribution being used.

   .. automethod:: solid_angle

      The default settings is to consider the entire solid angle. The optional
      *azimuth* and *elevation* arguments may be used to restrict the solid
      angle by specifying an interval of acceptable angular values, in deg.

   .. automethod:: target

      This method is best used in conjunction with a forward simulation, as it
      generates :ref:`target points <montecarlo:Target point>` over the surface
      of the specified *box*.

      The generation method is dependent on the direction configuration. If the
      direction is :py:func:`fixed <direction>` (the azimuth and elevation
      values referring to the box origin), then all generated particles are
      colinear. Otherwise, the particles are ingoing on the box, with a cosine
      angular distribution (corresponding to an isotropic flux entering the
      box).

      .. note::

         In the event that the :py:func:`solid angle <solid_angle>` of incoming
         particles is restricted, then a rejection sampling method is employed.


----

.. autoclass:: danton.Physics

   This class provides an interface to the Monte Carlo physics, which may be
   configured through attributes. See the :doc:`physics` section for further
   details.

   .. autoattribute:: bremsstrahlung

      The possible values of bremsstralung models are summarised in
      :numref:`tab-bremsstrahlung` below, the default setting being
      :python:`"SSR"`.

      .. _tab-bremsstrahlung:

      .. list-table:: Available bremsstrahlung models.
         :width: 75%
         :widths: auto
         :header-rows: 1

         * - Model
           - Reference
         * - :python:`"ABB"`
           - Andreev, Bezrukov and Bugaev, Physics of Atomic Nuclei 57 (1994)
             2066.
         * - :python:`"KKP"`
           - Kelner, Kokoulin and Petrukhin, Moscow Engineering Physics Inst.,
             Moscow, 1995.
         * - :python:`"SSR"`
           - [SSR19]_.

   .. autoattribute:: dis

      The possible values of predefined cross-section models are summarised in
      :numref:`tab-dis` below, the default setting being :python:`"CMS11"`.

      .. _tab-dis:

      .. list-table:: Available DIS models.
         :width: 75%
         :widths: auto
         :header-rows: 1

         * - Model
           - Reference
         * - :python:`"BGR18"`
           - [BGR18]_.
         * - :python:`"CMS11"`
           - [CMS11]_.
         * - :python:`"LO"`
           - Leading Order computation using the selected :py:attr:`PDFs <pdf>`
             (as in R. Gandhi et al., Astropart.Phys.5:81-110 (1996)).

      .. topic:: Custom cross-sections

         Alternatively to predefined models, one might specify a path to a text
         file containing cross-section values in *exactly* the following format
         (:download:`BGR18.txt <interface/BGR18.txt>`).

   .. autoattribute:: pair_production

      The possible values of pair-production models are summarised in
      :numref:`tab-pair-production` below, the default setting being
      :python:`"SSR"`.

      .. _tab-pair-production:

      .. list-table:: Available pair-production models.
         :width: 75%
         :widths: auto
         :header-rows: 1

         * - Model
           - Reference
         * - :python:`"KKP"`
           - Kelner, Kokoulin and Petrukhin, Soviet Journal of Nuclear Physics 7
             (1968) 237.
         * - :python:`"SSR"`
           - [SSR19]_.

   .. autoattribute:: pdf

   .. autoattribute:: photonuclear

      The possible values of photonuclear interaction models are summarised in
      :numref:`tab-photonuclear` below, the default setting being
      :python:`"DRSS"`.

      .. _tab-photonuclear:

      .. list-table:: Available photonuclear interaction models.
         :width: 75%
         :widths: auto
         :header-rows: 1

         * - Model
           - Reference
         * - :python:`"BBKS"`
           - Bezrukov, Bugaev, Sov. J. Nucl. Phys. 33 (1981), 635, with improved
             photon-nucleon cross-section according to `Kokoulin`_ and hard
             component from `Bugaev and Shlepin`_.
         * - :python:`"BM"`
           - Butkevich and Mikheyev, Soviet Journal of Experimental and
             Theoretical Physics 95 (2002) 11.
         * - :python:`"DRSS"`
           - Dutta, Reno, Sarcevic and Seckel, Phys.Rev. D63 (2001) 094020.

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
      :external:py:class:`ndarray <numpy.ndarray>`, each).

----

.. autofunction:: danton.compute

----

.. autofunction:: danton.particles

   This function returns a `structured <StructuredArray_>`_
   :external:py:class:`ndarray <numpy.ndarray>`  with the given *shape*.
   Particles are initialised with default properties, if not overriden by
   specifying *kwargs*. For instance, the following creates an array of 100
   particles (taus, by default) with an energy of 1 EeV, starting from the
   Equator (by default), and going to the North.

   >>> particles = danton.particles(100, energy=1E+09)

   The data structure (:external:py:class:`dtype <numpy.dtype>`) of a particle
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

.. _builder: https://en.wikipedia.org/wiki/Builder_pattern
.. _Bugaev and Shlepin: https://doi.org/10.1103/PhysRevD.67.034027
.. _ECEF: https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system
.. _EGM96: https://cddis.nasa.gov/926/egm96/egm96.html
.. _ellipsoid: https://en.wikipedia.org/wiki/Earth_ellipsoid
.. _geoid: https://en.wikipedia.org/wiki/Geoid
.. _Kokoulin: https://doi.org/10.1016/S0920-5632(98)00475-7
.. _PREM: https://en.wikipedia.org/wiki/Preliminary_reference_Earth_model
.. _PdgScheme: https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
.. _SRTMGL1: https://lpdaac.usgs.gov/products/srtmgl1v003/
.. _standard rock: https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/standard_rock.html
.. _StructuredArray: https://numpy.org/doc/stable/user/basics.rec.html
.. _TOML: https://toml.io/en/
.. _WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
