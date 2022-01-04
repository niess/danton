# Danton ![Build](https://github.com/niess/danton/workflows/Build/badge.svg)
( **D**ec**A**yi**N**g **T**aus fr**O**m **N**eutrinos )

## Description
Danton is a __C99 library__ dedicated to the sampling of decaying taus from
ultra high energy neutrinos interacting in the Earth. It can run in forward or
backward Monte Carlo. it can also be configured to sample tau fluxes instead of
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
Danton requires a UNIX system, eg. Linux or OSX. On Linux, the latest version of
the danton executable is available as an AppImage, from the
[releases](https://github.com/niess/danton/releases) assets, e.g. as:
```bash
# Fetch the danton AppImage from GitHub
wget https://github.com/niess/danton/releases/download/v0.3/danton-x86_64.AppImage

# Make the downloaded binary executable
chmod u+x danton-x86_64.AppImage

# Run the AppImage
./danton-x86_64.AppImage path/to/card.json
```

Alternatively, the provided [Makefile](Makefile) builds `libdanton` as a shared
library as well as the `danton` executable, e.g. as:
```bash
# Fetch the source and deps from GitHub
git clone --recursive https://github.com/niess/danton.git

# Build Danton
cd danton
make

# Run the executable
./bin/danton path/to/card.json
```

Note that `gfortran` is required in order to compile TAUOLA.

## API documentation
A documentation of the `libdanton` API is available [online][API:docs].

[API:docs]: https://niess.github.io/danton-docs

## Steering files

A syntaxic summary of the steering files parameters is provided here. Examples
can also be found in the [examples/cards](examples/cards) folder.

### Root items
```
decay           boolean              If `true` the sampled taus are decayed.
events          integer              The number of tentative Monte Carlo events.
longitudinal    boolean              If `true` the transverse transport is disabled.
mode            string               The run mode, one of "backward", "forward" or "grammage".
output-file     string, null         The output file name or `null` for `stdout`.
requested       integer              The requested number of valid Monte Carlo events.
seed            unsigned, null       The simulation random seed. If null, then the seed is set
                                       from the OS entropy using /dev/urandom.
```

In addition to the previous general parameters one also has the following keys :
`"earth-model"`, `"particle-sampler"`, `"physics"`, `"primary-flux"`, and
`"stepping"`. The corresponding options are described hereafter.

### Earth model
```
reference       string               The reference model: "PREM" (spherical),
                                       "WGS84" (elliptical) or "EGM96" (geoidal).
sea             boolean              If `true` the Earth is covered with sea.
topography      string, float        Path to a folder containing topography data, e.g.
                                       SRTMGL1.v3 tiles. Alternatively a float value can be
                                       provided, specifying a constant topography altitude.
material        string               The material composing the topography.
density         float                The density of the topography, in kg/m^3.
```

Note that the legacy PREM has an external layer of 3km of sea water with a
density of 1.02 g/cm</sup>3</sup>. If the sea is disabled this layer is replaced
with [Standard Rock][1].

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

### Physics
```
bremsstrahlung  string, null         Model for the bremsstrahlung process, for taus.
pair-production string, null         Model for the pair production process, for taus.
photonuclear    string, null         Model for photonuclear interactions, for taus.
DIS             string, null         Model for Deep Inelastic Scattering (DIS),
                                       for neutrinos.
PDF             string, null         Use a specific PDF set for neutrinos DIS.
```

For tau energy losses (`"bremsstrahlung"`, `"pair-production"` and
`"photonuclear"`) the model must correspond to one available in
[PUMAS](https://pumas.readthedocs.io/en/latest/api/#HEAD/type/pumas_physics_settings).

For DIS cross-sections, the supported models are `"LO"` (Leading Order
computation, see e.g. [Gandhi et al.,
(1995)](https://arxiv.org/abs/hep-ph/9512364)), `"CSMS"` ([Cooper-Sarkar,
Mertsch and Sarkar](https://arxiv.org/abs/1106.3723)).  and  `"BGR18"`
([Bertone, Gold and Rojo (2018)](https://arxiv.org/abs/1808.02034)).
Alternatively, one can provide a path to a file containing cross-section values
in [ENT](https://github.com/niess/ent)'s format.

For PDF sets, the built-in models are `"HERAPDF15NLO"`, `"CT14nlo"` and
`"NNPDF31sx"`.  Alternatively, one can provide a path to a file containing
cross-section values in Les Houches Accord (LHA) format.


_Note that `"physics"` options must be set **before** `"earth"` ones in the
data-card._

### Primary flux
```
$particle       [$model, {...}]      The primary spectrum model for the corresponding
                                       particle.
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
append          boolean              If `true`, append to an existing stepping file,
                                       instead of overwriting.
path            string               Path to the stepping file.
```

The stepping file contains a summary of Monte Carlo steps, in JSON format.

## License
The Danton source is distributed under the **GNU LGPLv3** license. See the
provided [LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
