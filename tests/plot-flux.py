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
    """Plot the Monte Carlo flux results."""

    data = []
    for path in args.path:
        with gzip.open(path, "rb") as f:
            data.append(pickle.load(f))

    plt.figure()
    for i, path in enumerate(args.path):
        clr = "k" if "backward" in path else "r"
        energy = Histogram.new(data[i], "energy")
        energy.errorbar(fmt=f"{clr}o")
    plt.xlabel("energy (GeV)")
    plt.ylabel("flux (GeV$^{-1}$)")
    plt.ylim(1E-23, 1E-11)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a flux simulation")
    parser.add_argument("path",
        help = "Input simulation results",
        nargs = "+"
    )

    args = parser.parse_args()
    plot(args)
