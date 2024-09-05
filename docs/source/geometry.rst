Geometry description
====================

Earth structure
---------------

Danton Model's the Earth interior structure according to the Preliminary
Reference Earth Model (`PREM`_). That is, the Earth is composed of concentric
layers, ranging from the very core to the surface, with varying bulk densities
depending on the depth.

.. note::

   The legacy `PREM`_ Earth is covered by a 3km deep ocean, without any land.
   This can be disabled however (see the :py:attr:`ocean
   <danton.Geometry.ocean>` attribute of the :py:class:`Geometry
   <danton.Geometry>` class).

The PREM originates from seismic data. It does not specify the Earth's inner
composition. Currently, Danton assumes a uniform composition corresponding to
`standard rock`_, except for the optional ocean layer.


Topography
----------

Atmosphere
----------


.. ============================================================================
.. 
.. URL links.
.. 
.. ============================================================================

.. _PREM: https://en.wikipedia.org/wiki/Preliminary_reference_Earth_model
.. _standard rock: https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/standard_rock.html
