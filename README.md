# DANTON
( **D**ec**A**yi**N**g **T**aus fr**O**m **N**eutrinos )

## Description
DANTON is a __C99 library__ dedicated to the sampling of decaying taus from
ultra high energy neutrinos interacting in the Earth. It can run in forward or
backward Monte-Carlo. it can also be configured to sample tau fluxes instead of
decay densities, or to sample transmitted neutrinos fluxes.

The library is shipped with an __executable__, `danton` which takes a
*data card* in __JSON__ format as argument, e.g. :
```bash
danton mycard.json
```

Currently a spherical Earth is assumed with a density given by the
**P**reliminary **E**arth **M**odel (**PEM**). The atmosphere is
**US standard**.

## Installation
Currently there is no automatic build procedure. On a linux box you might try
the provided [setup.sh](setup.sh) and [Makefile](Makefile) as :
```bash
. setup.sh
make lib && make
```
This will build the dynamic libraries for all the submodules and then the
`libdanton.so` library and the `danton` executable.

## Data cards

A syntaxic summary of the data cards options is provided here. Examples can
also be found in the [cards](cards) folder.

### Root keys
```
decay           boolean              If `true` the sampled taus are decayed.
events          integer              The number of Monte-Carlo events to run.
longitudinal    boolean              If `true` the transverse transport is disabled.
mode            string               The run mode, one of "backward", "forward" or "grammage".
output-file     string, null         The output file name or `null` for `stdout`.
```

In addition to the previous general options one also has the three following
keys : `"earth-model"`, `"particle-sampler"` and `"primary-flux"`. The
corresponding options are described hereafter.

### Earth model
```
sea             boolean              If `true` the PEM Earth is covered with sea.
```

Note that the legacy PEM has an external layer of 3km of sea water. If the sea
is disabled this layer is replaced with [Standard Rock][1].

[1]: http://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/standard_rock.html

### Particle sampler
```
altitude        float, float[2]      The altitude (range) of the sampled particles.
elevation       float, float[2]      The elevation angle (range) of the sampled particles.
energy          float, float[2]      The energy (range) of the sampled particles.
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

## License
The sources are under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files.
