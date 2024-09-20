#! /usr/bin/env python3
import argparse
import danton
import gzip
import numpy
from pathlib import Path
import pickle


PREFIX = Path(__file__).parent

import sys
sys.path.append(str(PREFIX))
from histogram import Histogram


def run(args):
    """Run an effective area computation."""

    simulation = danton.Simulation(
        mode = "backward",
        tau_decays = True,
        longitudinal = True,
        geoid = "WGS84",
        ocean = args.ocean
    )

    if args.topography is not None:
        simulation.geometry.topography = args.topography

    n = args.number_of_events
    emin, emax = args.energy_min, args.energy_max

    box = simulation.box(
        size = [args.width, args.length, args.height],
        latitude = args.latitude,
        longitude = args.longitude,
        altitude = args.altitude
    )

    generator = simulation.particles()      \
        .pid(15)                            \
        .inside(box, limit=args.limit)      \
        .powerlaw(emin, emax, exponent=-1)

    if args.direction is not None:
        assert(args.azimuth is None)
        assert(args.elevation is None)

        generator.direction(args.direction)

    elif (args.azimuth is not None) or (args.elevation is not None):
        assert(args.direction is None)
        assert(args.azimuth is not None)
        assert(args.elevation is not None)

        generator.direction(
            azimuth = args.azimuth,
            elevation = args.elevation
        )

    else:
        generator.solid_angle(elevation=[-10, 10])

    secondaries, N = generator.generate(n)

    result = simulation.run(secondaries)

    primaries = result.primaries
    sel = (primaries["pid"] == 16) & (primaries["energy"] <= emax)
    primaries = primaries[sel]

    data = {
        "seed": simulation.random.seed,
        "version": danton.VERSION
    }
    data.update(args.__dict__)

    if args.all:
        secondaries = secondaries[primaries["event"]]
        data["events"] = N
        data["primaries"] = primaries
        data["secondaries"] = secondaries
    else:
        tmp = {
            "events": N, "energy_min": emin, "energy_max": emax,
            "primaries": primaries
        }
        data["area"] = Histogram.new(tmp, "energy", particles="primaries")

    outfile = f"area-{simulation.random.seed:0X}.pkl.gz"
    if args.output_directory is None:
        path = PREFIX / "data"
    else:
        path = Path(args.output_directory)
    path.mkdir(parents=True, exist_ok=True)
    outfile = path / outfile

    with gzip.open(outfile, "wb") as f:
        pickle.dump(data, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Run an effective area computation"
    )
    parser.add_argument("-n", "--number-of-events",
        help = "Number of tentative secondaries",
        type = int,
        default = 10000
    )
    parser.add_argument("--limit",
        help = "Use the distance limit preselection",
        action = argparse.BooleanOptionalAction,
        default = True
    )

    position = parser.add_argument_group("Position",
        description = "Options controlling the box position."
    )
    position.add_argument("--latitude",
        help = "Box center's latitude, in deg",
        type = float,
        default = 0
    )
    position.add_argument("--longitude",
        help = "Box center's longitude, in deg",
        type = float,
        default = 0
    )
    position.add_argument("--altitude",
        help = "Box center's altitude, in m",
        type = float,
        default = 0
    )

    direction = parser.add_argument_group("Direction",
        description = "Options controlling the direction."
    )
    direction.add_argument("--azimuth",
        help = "Azimuth angle of the direction, in deg",
        type = float,
    )
    direction.add_argument("--elevation",
        help = "Elevation angle of the direction, in deg",
        type = float,
    )
    direction.add_argument("--direction",
        help = "Geocentric direction (cartesian coordinates)",
        type = float,
        nargs = 3,
        metavar = ("X", "Y", "Z"),
    )

    size = parser.add_argument_group("Size",
        description = "Options controlling the box size."
    )
    size.add_argument("--width",
        help = "Box width along the x-axis, in m",
        type = float,
        default = 1E+05
    )
    size.add_argument("--length",
        help = "Box length along the y-axis, in m",
        type = float,
        default = 1E+05
    )
    size.add_argument("--height",
        help = "Box height along the z-axis, in m",
        type = float,
        default = 1E+04
    )

    energy = parser.add_argument_group("Energy",
        description = "Options controlling the sampled energy range."
    )
    energy.add_argument("--energy-min",
        help = "Minimum energy, in GeV",
        type = float,
        default = 1E+06
    )
    energy.add_argument("--energy-max",
        help = "Maximum energy, in GeV",
        type = float,
        default = 1E+12
    )

    geometry = parser.add_argument_group("Geometry",
        description = "Options controlling the Earth geometry."
    )
    geometry.add_argument("--ocean",
        help = "Enable PREM ocean (when no topography)",
        action = "store_true",
        default = False
    )
    geometry.add_argument("-t", "--topography",
        help = "Path to topography data"
    )

    output = parser.add_argument_group("Output",
        description = "Options controlling the output data."
    )
    output.add_argument("-a", "--all",
        help = "Log Monte Carlo events, instead of a summary histogram.",
        action = "store_true"
    )
    output.add_argument("-o", "--output-directory",
        help = "Output directory"
    )

    args = parser.parse_args()
    run(args)
