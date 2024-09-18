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
    """Plot an effective area estimate."""

    data = []
    for i, path in enumerate(args.path):
        with gzip.open(path, "rb") as f:
            d = pickle.load(f)

            if i == 0:
                if (d["azimuth"] is None) and \
                   (d["elevation"] is None) and \
                   (d["direction"] is None):
                    unit = f"{args.unit} sr"
                else:
                    unit = args.unit

            data.append(d["area"])

    if args.unit == "cm":
        scale = 1E+04
    elif args.unit == "m":
        scale = 1
    elif args.unit == "km":
        scale = 1E-06
    else:
        raise NotImplementedError()

    area = Histogram \
        .sum(data) \
        .scaled(yscale=scale)

    plt.figure()
    area.errorbar(fmt=f"ko")
    plt.xlabel("energy (GeV)")
    plt.ylabel(f"area ({unit}$^2$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot an effective area estimate")
    parser.add_argument("path",
        help = "Input simulation results",
        nargs = "+"
    )
    parser.add_argument("-u", "--unit",
        help = "Unit for the effective area plot",
        choices = ["cm", "m", "km"],
        default = "m"
    )

    args = parser.parse_args()
    plot(args)
