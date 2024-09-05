Danton
======
*(DecAyiNg Taus frOm Neutrinos)*

----

Danton is a Python package designed specifically for the sampling of tau decays
from ultra-high energy neutrinos interacting with the Earth's surface. It is
capable of running in forward or backward Monte Carlo, and can also be
configured to sample tau fluxes instead of decay densities, or to sample
transmitted neutrino fluxes.

The interface has been designed with simplicity in mind. That is, Monte Carlo
:py:func:`particles <danton.particles>` are :py:meth:`injected
<danton.Simulation.run>` into the simulation :py:class:`Geometry
<danton.Geometry>` as a :external:py:class:`numpy.ndarray`, and
:external:py:class:`numpy.ndarray`\ s of Monte Carlo states are returned. This
basic workflow is illustrated below,

.. code:: python

   import danton

   simulation = danton.Simulation()
   particles = danton.particles(10000, energy=1E+09, elevation=1.0)
   result = simulation.run(particles)


System of units
---------------

.. note::

   Danton employs the Metre-Kilogram-Second (MKS) system of units (e.g. kg/m\
   :sup:`3` for a density), with the exception of energies and momenta, which
   are expressed in GeV and GeV/c, respectively. Additionally, angles (e.g.
   geodetic coordinates) are expressed in degrees.

Documentation
-------------

.. toctree::
   :maxdepth: 2

   geometry
   materials
   interface
