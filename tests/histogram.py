import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple


class Histogram(NamedTuple):
    """Histogramed data."""

    x: np.ndarray
    y: np.ndarray
    xerr: np.ndarray
    yerr: np.ndarray

    @classmethod
    def new(cls, data, attribute):
        """Create a new histogram from Monte Carlo data."""

        n = data["events"]
        secondaries = data["secondaries"]
        samples, weights = secondaries[attribute], secondaries["weight"]

        if attribute == "energy":
            edges = np.geomspace(data["energy_min"], data["energy_max"], 61)
            widths = edges[1:] - edges[:-1]
            x = np.sqrt(edges[1:] * edges[:-1])
            xerr = np.array([x - edges[:-1], edges[1:] - x])
        else:
            if attribute == "elevation":
                edges = np.linspace(
                    -data["elevation"],
                    data["elevation"],
                    41
                )
            elif attribute == "altitude":
                edges = np.linspace(
                    0,
                    data["height"] * 3 / 4,
                    41
                )
            else:
                raise NotImplementedError(attribute)

            widths = edges[1:] - edges[:-1]
            x = 0.5 * (edges[1:] + edges[:-1])
            xerr = np.array([x - edges[:-1], edges[1:] - x])

        y, _ = np.histogram(samples, edges, weights=weights)
        yerr, _ = np.histogram(samples, edges, weights=weights**2)
        y /= (n * widths)
        yerr = np.sqrt((yerr / n - y**2) / n) / widths

        return cls(x, y, xerr, yerr)


    def errorbar(self, **kwargs):
        """Draw the histogram using errorbars."""

        plt.errorbar(self.x, self.y, xerr=self.xerr, yerr=self.yerr, **kwargs)

    def plot(self, *args, **kwargs):
        """Draw a plot of the histogram."""

        plt.plot(self.x, self.y, *args, **kwargs)
