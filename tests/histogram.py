import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple


class Histogram(NamedTuple):
    """Histogramed data."""

    x: np.ndarray
    y: np.ndarray
    xerr: np.ndarray
    yerr: np.ndarray
    n: int

    @classmethod
    def new(cls, data, attribute, particles="secondaries", raw=False):
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
        if not raw:
            y /= (n * widths)
            yerr = np.sqrt(np.maximum(yerr / n - y**2, 0.0) / n) / widths

        return cls(x, y, xerr, yerr, n)


    @classmethod
    def sum(cls, histograms, raw=False):
        """Sum of histograms."""

        x = histograms[0].x
        for h in histograms[1:]:
            assert((h.x == x).all())
        xerr = histograms[0].xerr

        if raw:
            y, yerr, n = 0, 0, 0
            for h in histograms:
                y += h.y
                yerr += h.yerr
                n += h.n

            widths = xerr[1,:] + xerr[0,:]
            y /= (n * widths)
            yerr = np.sqrt(np.maximum(yerr / n - y**2, 0.0) / n) / widths

        else:
            y = np.empty(x.size)
            yerr = np.empty(x.size)
            for i in range(x.size):
                w = np.array([
                    1 / h.yerr[i]**2 if h.yerr[i] > 0.0 else 0.0
                        for h in histograms
                ])
                yi = np.array([h.y[i] for h in histograms])
                tmp = sum(w)
                if tmp > 0:
                    y[i] = sum(yi * w) / tmp
                    yerr[i] = 1 / np.sqrt(tmp)
                else:
                    y[i] = 0.0
                    yerr[i] = 0.0

            n = sum(h.n for h in histograms)

        return cls(x, y, xerr, yerr, n)


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

        return self.__class__(x, y, xerr, yerr, self.n)

    def to_raw(self):
        """Convert to raw representation."""

        x, xerr, n = self.x, self.xerr, self.n
        widths = xerr[1,:] + xerr[0,:]
        y = self.y * n * widths
        yerr = ((self.yerr * widths)**2 * n + self.y**2) * n

        return self.__class__(x, y, xerr, yerr, self.n)
