#! /usr/bin/env python3
import argparse
import danton
import gzip
import numpy
from pathlib import Path
import pickle
import psutil


PREFIX = Path(__file__).parent


def run(args):
    """Run a flux simulation."""

    simulation = danton.Simulation(
        mode = args.mode,
        tau_decays = True,
        longitudinal = True,
        geoid = "PREM"
    )

    n = args.number_of_events
    emin, emax = args.energy_min, args.energy_max

    def primary_flux(e):
        return (emin * emax) / (e**2 *(emax - emin))

    box = simulation.box(
        size = [args.width, args.width, args.height],
        altitude = args.height / 4
    )

    process = psutil.Process()

    def get_cpu():
        t = process.cpu_times()
        return t.user + t.system

    cpu = get_cpu()

    if args.mode == "backward":
        particles, size = simulation.particles(weight=True)           \
            .pid(15)                                                  \
            .powerlaw(emin, emax, exponent=-1)                        \
            .inside(box)                                              \
            .solid_angle(elevation=[-args.elevation, args.elevation]) \
            .generate(n)
    else:
        particles, size = simulation.particles(weight=True)           \
            .pid(16)                                                  \
            .powerlaw(emin, emax, exponent=-1)                        \
            .target(box)                                              \
            .solid_angle(elevation=[-args.elevation, args.elevation]) \
            .generate(n)

    result = simulation.run(particles)
    cpu = get_cpu() - cpu

    if args.mode == "backward":
        primaries = result.primaries
        sel = (primaries["pid"] == 16) & (primaries["energy"] <= emax)
        primaries = primaries[sel]
        secondaries = particles[primaries["event"]]
        secondaries["weight"] = primary_flux(primaries["energy"]) * \
            primaries["weight"]
        random_index = primaries["random_index"]
    else:
        secondaries = result.secondaries
        sel = (secondaries["pid"] == 15) & (secondaries["energy"] >= emin) & \
            (secondaries["elevation"] >= -args.elevation) & \
            (secondaries["elevation"] <= args.elevation) & \
            box.inside(secondaries)
        secondaries = secondaries[sel]
        primaries = particles[secondaries["event"]]
        secondaries["weight"] *= primary_flux(primaries["energy"])
        random_index = secondaries["random_index"]

    data = {
        "mode": args.mode,
        "events": size,
        "elevation": args.elevation,
        "width": args.width,
        "height": args.height,
        "energy_min": emin,
        "energy_max": emax,
        "secondaries": secondaries,
        "seed": simulation.random.seed,
        "random_index": random_index,
        "cpu": cpu,
        "version": danton.VERSION
    }

    tag = "_".join([
        f"{args.mode}",
        f"{args.elevation:.3f}",
        f"{args.width:.0E}",
        f"{args.height:.0E}",
        f"{simulation.random.seed:0X}"
    ])
    outfile = f"decays-{tag}.pkl.gz"
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
    parser.add_argument("-H", "--height",
        help = "Box height, in m",
        type = float,
        default = 2E+02
    )
    parser.add_argument("-W", "--width",
        help = "Box width, in m",
        type = float,
        default = 1E+04
    )
    parser.add_argument("-e", "--elevation",
        help = "Elevation half-width, in deg",
        type = float,
        default = 5
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
