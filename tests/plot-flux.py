#! /usr/bin/env python3
import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pickle


PREFIX = Path(__file__).parent

plt.style.use(PREFIX / "paper.mplstyle")


def histogram(path):
    """Histogram flux data."""

    with gzip.open(path, "rb") as f:
        data = pickle.load(f)

    n = data["events"]
    secondaries = data["secondaries"]
    e, w = secondaries["energy"], secondaries["weight"]

    edges = np.geomspace(data["energy_min"], data["energy_max"], 61)
    widths = edges[1:] - edges[:-1]
    x = np.sqrt(edges[1:] * edges[:-1])
    xerr = [x - edges[:-1], edges[1:] - x]

    y, _ = np.histogram(e, edges, weights=w)
    yerr, _ = np.histogram(e, edges, weights=w**2)
    y /= (n * widths)
    yerr = np.sqrt((yerr / n - y**2) / n) / widths

    return x, y, xerr, yerr


def plot(args):
    """Plot the Monte Carlo flux results."""

    plt.figure()
    for path in args.path:
        clr = "k" if "backward" in path else "r"
        x, y, xerr, yerr = histogram(path)
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=f"{clr}o")

    plt.xlabel("energy (GeV)")
    plt.ylabel("flux (GeV$^-1$)")
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
