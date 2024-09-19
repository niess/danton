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

    if args.unit == "cm":
        scale = 1E+04
    elif args.unit == "m":
        scale = 1
    elif args.unit == "km":
        scale = 1E-06
    else:
        raise NotImplementedError()

    def plot_data(paths, fmt=None, label=None):
        data = []
        for i, path in enumerate(paths):
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

        area = Histogram \
            .sum(data) \
            .scaled(yscale=scale)

        if fmt is None:
            fmt = "ko"

        area.errorbar(fmt=fmt, label=label)
        return unit

    plt.figure()
    legend = False
    if args.path:
        unit = plot_data(args.path)
    else:
        for data in args.label:
            label, fmt, *paths = data
            unit = plot_data(paths, fmt, label)
            legend = True

    plt.xlabel("energy (GeV)")
    plt.ylabel(f"effective area ({unit}$^2$)")
    plt.xscale("log")
    plt.yscale("log")
    if args.ylim:
        plt.ylim(*args.ylim)
    if legend:
        plt.legend()
    if args.save_fig:
        plt.savefig(args.save_fig)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot an effective area estimate")
    parser.add_argument("path",
        help = "Input simulation results",
        nargs = "*"
    )
    parser.add_argument("-l", "--label",
        help = "Labelled simulation results",
        nargs = "+",
        action = "append"
    )
    parser.add_argument("-s", "--save-fig",
        help = "Save the figure as PATH",
        metavar = "PATH"
    )

    scale = parser.add_argument_group("Scale",
        "Options controlling the plot scale."
    )
    scale.add_argument("-u", "--unit",
        help = "Unit for the effective area plot",
        choices = ["cm", "m", "km"],
        default = "m"
    )
    scale.add_argument("-y", "--ylim",
        help = "Limits along the y-axis",
        type = float,
        nargs = 2
    )

    args = parser.parse_args()
    plot(args)
