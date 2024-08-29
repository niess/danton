#! /usr/bin/env python3
import danton
import numpy
import matplotlib.pyplot as plot


# Instanciate the Monte Carlo geometry.
geometry = danton.Geometry(geoid="EGM96")

# Prepare the coordinates grid.
latitude = numpy.linspace(-90, 90, 181)
longitude = numpy.linspace(-180, 180, 361)
lon, lat = numpy.meshgrid(longitude, latitude)

# Get geoid undulations.
undulations = geometry.geoid_undulation(latitude=lat, longitude=lon)

# Plot the result.
plot.figure()
plot.pcolormesh(
    longitude,
    latitude,
    undulations,
    cmap = "seismic",
    vmin = -100,
    vmax = 100
)
plot.colorbar()
plot.xlabel("longitude (deg)")
plot.ylabel("latitude (deg)")
plot.title("Geoid undulation (m)")
plot.show()
