#! /usr/bin/env python3
import argparse
import danton
import gzip
import numpy
from pathlib import Path
import pickle


PREFIX = Path(__file__).parent


def run(args):
    """Run a flux simulation."""

    simulation = danton.Simulation(
        mode = args.mode,
        tau_decays = False,
        longitudinal = True,
        geoid = "PREM"
    )

    n = args.number_of_events
    emin, emax = args.energy_min, args.energy_max

    def primary_flux(e):
        return (emin * emax) / (e**2 *(emax - emin))

    pid = 15 if args.mode == "backward" else 16

    particles = simulation.particles(weight=True) \
        .pid(pid)                                 \
        .powerlaw(emin, emax, exponent=-1)        \
        .position(altitude=args.altitude)         \
        .direction(elevation=args.elevation)      \
        .generate(n)

    result = simulation.run(particles)

    if args.mode == "backward":
        primaries = result.primaries
        sel = (primaries["pid"] == 16) & (primaries["energy"] <= emax)
        primaries = primaries[sel]
        secondaries = particles[primaries["event"]]
        secondaries["weight"] = primary_flux(primaries["energy"]) * \
            primaries["weight"]
    else:
        secondaries = result.secondaries
        sel = (secondaries["pid"] == 15) & (secondaries["energy"] >= emin)
        secondaries = secondaries[sel]
        primaries = particles[secondaries["event"]]
        secondaries["weight"] *= primary_flux(primaries["energy"])

    w = secondaries["weight"]
    mu = sum(w) / n
    sigma = numpy.sqrt((sum(w**2) / n - mu**2) / n)
    print(f"{mu:.3E} +- {sigma:.3E}")

    data = {
        "mode": args.mode,
        "events": n,
        "elevation": args.elevation,
        "altitude": args.altitude,
        "energy_min": emin,
        "energy_max": emax,
        "secondaries": secondaries
    }

    outfile = f"flux-{args.mode}-{args.elevation:.3f}-{args.altitude:.0E}.pkl.gz"
    if args.output_directory is None:
        path = PREFIX / "data"
    else:
        path = Path(args.output_directory)
    path.mkdir(parents=True, exist_ok=True)
    outfile = path / outfile

    with gzip.open(outfile, "wb") as f:
        pickle.dump(data, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a flux simulation")
    parser.add_argument("-m", "--mode",
        help = "Monte Carlo mode",
        choices = ["backward", "forward"],
        default = "backward"
    )
    parser.add_argument("-n", "--number-of-events",
        help = "Monte Carlo statistics",
        type = int,
        default = 10000
    )
    parser.add_argument("-a", "--altitude",
        help = "Sampling altitude, in m",
        type = float,
        default = 0
    )
    parser.add_argument("-e", "--elevation",
        help = "Sampling elevation angle, in deg",
        type = float,
        default = 1
    )
    parser.add_argument("--energy_min",
        help = "Minimum energy, in GeV",
        type = float,
        default = 1E+06
    )
    parser.add_argument("--energy_max",
        help = "Maximum energy, in GeV",
        type = float,
        default = 1E+12
    )
    parser.add_argument("-o", "--output-directory",
        help = "Output directory"
    )

    args = parser.parse_args()
    run(args)
