#! /usr/bin/env python3
import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle


PREFIX = Path(__file__).parent

import sys
sys.path.append(str(PREFIX))
from histogram import Histogram

plt.style.use(PREFIX / "paper.mplstyle")


def plot(args):
    """Plot the Monte Carlo decays results."""

    data = []
    for path in args.path:
        with gzip.open(path, "rb") as f:
            data.append(pickle.load(f))

    def draw(path, histogram):
        """Draw an histogram."""

        backward = "backward" in path
        label = "backward" if backward else "forward"
        clr = "k" if backward else "r"
        if args.mode == "plot" or (args.mode == "mixed" and backward):
            histogram.plot(f"{clr}-", label=label)
        else:
            histogram.errorbar(fmt=f"{clr}o", label=label)

    plt.figure()
    for i, path in enumerate(args.path):
        energy = Histogram.new(data[i], "energy")
        draw(path, energy)
    plt.xlabel("energy (GeV)")
    plt.ylabel("rate (GeV$^{-1}$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()

    plt.figure()
    for i, path in enumerate(args.path):
        elevation = Histogram.new(data[i], "elevation")
        draw(path, elevation)
    plt.xlabel("elevation (deg)")
    plt.ylabel("rate (deg$^{-1}$)")
    plt.yscale("log")
    plt.legend()

    plt.figure()
    for i, path in enumerate(args.path):
        altitude = Histogram.new(data[i], "altitude")
        draw(path, altitude)
    plt.xlabel("altitude (m)")
    plt.ylabel("rate (m$^{-1}$)")
    plt.yscale("log")
    plt.legend()

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a flux simulation")
    parser.add_argument("path",
        help = "Input simulation results",
        nargs = "+"
    )
    parser.add_argument("-m", "--mode",
        help = "Plot mode",
        choices = ["errorbar", "mixed", "plot"],
        default = "errorbar"
    )

    args = parser.parse_args()
    plot(args)
