# Danton [![Documentation Status](https://readthedocs.org/projects/danton/badge/?version=latest)](https://danton.readthedocs.io/en/latest/?badge=latest)
( **D**ec**A**yi**N**g **T**aus fr**O**m **N**eutrinos )


## Description

Danton is a Python package designed specifically for the sampling of tau decays
from ultra-high energy neutrinos interacting with the Earth's surface. It is
capable of running in forward or backward Monte Carlo, and can also be
configured to sample tau fluxes instead of decay densities, or to sample
transmitted neutrino fluxes.

The interface has been designed with simplicity in mind. That is, Monte Carlo
particles are injected into the simulation geometry as a `numpy.ndarray`, and
`numpy.ndarray`s of Monte Carlo states are returned. This basic workflow is
illustrated below,

```python
import danton

simulation = danton.Simulation()
particles = danton.particles(10000, energy=1E+09, elevation=1.0)
result = simulation.run(particles)
```


## License
The Danton source is distributed under the **GNU LGPLv3** license. See the
provided [LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
