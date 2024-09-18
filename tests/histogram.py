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
    def new(cls, data, attribute, particles="secondaries"):
        """Create a new histogram from Monte Carlo data."""

        n = data["events"]
        particles = data[particles]
        samples, weights = particles[attribute], particles["weight"]

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


    @classmethod
    def sum(cls, histograms):
        """Weighted sum of histograms."""

        x = histograms[0].x
        for h in histograms[1:]:
            assert((h.x == x).all())
        xerr = histograms[0].xerr

        y = np.empty(x.size)
        yerr = np.empty(x.size)
        for i in range(x.size):
            w = np.array([
                1 / h.yerr[i]**2 if h.yerr[i] > 0.0 else 0.0 for h in histograms
            ])
            yi = np.array([h.y[i] for h in histograms])
            tmp = sum(w)
            if tmp > 0:
                y[i] = sum(yi * w) / tmp
                yerr[i] = 1 / np.sqrt(tmp)
            else:
                y[i] = 0.0
                yerr[i] = 0.0

        return cls(x, y, xerr, yerr)


    def errorbar(self, **kwargs):
        """Draw the histogram using errorbars."""

        plt.errorbar(self.x, self.y, xerr=self.xerr, yerr=self.yerr, **kwargs)


    def plot(self, *args, **kwargs):
        """Draw a plot of the histogram."""

        plt.plot(self.x, self.y, *args, **kwargs)

    def scaled(self, xscale=None, yscale=None):
        """Generate a scaled histogram."""

        if xscale is None:
            xscale = 1
        if yscale is None:
            yscale = 1

        yscale /= xscale

        x = self.x * xscale
        y = self.y * yscale
        xerr = self.xerr * xscale
        yerr = self.yerr * yscale

        return self.__class__(x, y, xerr, yerr)
