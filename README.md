[![Build Status](https://travis-ci.com/niess/danton.svg?branch=master)](https://travis-ci.com/niess/danton)

# DANTON
( **D**ec**A**yi**N**g **T**aus fr**O**m **N**eutrinos )

## Description
DANTON is a __C99 library__ dedicated to the sampling of decaying taus from
ultra high energy neutrinos interacting in the Earth. It can run in forward or
backward Monte-Carlo. it can also be configured to sample tau fluxes instead of
decay densities, or to sample transmitted neutrinos fluxes.

The library is shipped with an __executable__, `danton` which takes a *data
card* in __JSON__ format as argument, e.g. :
```bash
danton mycard.json
```

The Earth core is modelled according to the **P**reliminary **R**eference
**E**arth **M**odel ([**PREM**][PREM]). The atmosphere is **US standard**. A
detailed topography can be provided from world wide elevation models, e.g.
[ASTER][ASTER] or [SRTMGL1][SRTMGL1].

[PREM]: https://www.sciencedirect.com/science/article/pii/0031920181900467
[ASTER]: https://lpdaac.usgs.gov/dataset_discovery/aster/aster_products_table/astgtm_v002
[SRTMGL1]: https://lpdaac.usgs.gov/dataset_discovery/measures/measures_products_table/srtmgl1_v003

## Installation
DANTON has been tested on both Linux and OSX. Not on Windows though. The
provided [Makefile](Makefile) builds `libdanton` as a shared library as well
as the `danton` executable, e.g. as:
```bash
# Fetch the source and deps from GitHub
git clone --recursive https://github.com/niess/danton.git

# Build DANTON
cd danton
make

# Run the executable
./bin/danton path/to/card.json
```

DANTON requires `gfortran` in order to compile TAUOLA. In addition, the default Makefile
expects `libpng` and `libtiff` to be installed. Those are required by some topography
models. They can be disabled by editing the `USE_PNG` and `USE_TIFF` flags in the
Makefile.
 
## API documentation
A documentation of the `libdanton` API is available [online][API:docs].

[API:docs]: https://niess.github.io/danton-docs

## Steering files

A syntaxic summary of the steering files parameters is provided here. Examples
can also be found in the [share/cards](share/cards) folder.

### Root items
```
decay           boolean              If `true` the sampled taus are decayed.
events          integer              The number of Monte-Carlo events to run.
longitudinal    boolean              If `true` the transverse transport is disabled.
mode            string               The run mode, one of "backward", "forward" or "grammage".
output-file     string, null         The output file name or `null` for `stdout`.
requested       integer              The requested number of valid Monte-Carlo events
```

In addition to the previous general parameters one also has the following keys :
`"earth-model"`, `"particle-sampler"`, `"primary-flux"` and `"stepping"`. The
corresponding options are described hereafter.

### Earth model
```
geodesic        string               The geodesic model: "PREM" (spherical) or "WGS84".
sea             boolean              If `true` the PREM Earth is covered with sea.
topography      [string, integer]    The topography data location and the in-memory stack size.
```

Note that the legacy PREM has an external layer of 3km of sea water. If the sea
is disabled this layer is replaced with [Standard Rock][1].

[1]: http://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/standard_rock.html

### Particle sampler
```
altitude        float, float[2]      The altitude (range) of the sampled particles.
azimuth         float, float[2]      The azimuth angle (range) of the sampled particles.
elevation       float, float[2]      The elevation angle (range) of the sampled particles.
energy          float, float[2]      The energy (range) of the sampled particles.
latitude        float                The geodetic latitude of the sampled particles.
longitude       float                The geodetic longitude of the sampled particles.
weight          {$particle:float}    The name and weight of the particles to sample.
```

The valid *particle* names are `"nu_tau"`, `"nu_tau~"`, `"nu_mu"`, `"nu_mu~"`,
`"nu_e"`, `"nu_e~"`, `"tau"` and `"tau~"`.

### Primary flux
```
$particle       [$model, {...}]      The primary spectrum model for the corresponding particle.
```

Note that for a primary flux `$particle` can't be a `"tau"` or `"tau~"`. The
valid primary models are `"discrete"` and `"power-law"`, as described hereafter.

#### Discrete spectrum
```
energy          float                The total energy of the primary.
weight          float                The weight of the primary, i.e. the integrated flux.
```

#### Power law spectrum
```
energy          float[2]             The energy range of the primary spectrum.
exponent        float                The exponent of the power law.
weight          float                The weight of the primary, i.e. the integrated flux.
```

### Stepping
```
append          boolean              If `true`, append to the output file.
path            string               Path to the output file.
verbosity       integer              Verbosity level for recording Monte-Carlo steps.
```

## License
The sources are under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
