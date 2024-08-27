#! /usr/bin/env python3
import danton

# Create a Danton simulation interface.
simulation = danton.Simulation()

# Configure the Monte Carlo for backward sampling of tau decays.
simulation.mode = "backward"
simulation.tau_decays = True

# Set an elliptical Earth (using the WGS84 ellipsoid).
simulation.geometry.geoid = "WGS84"

# Generate a bunch of taus at decay, with an elevation angle of 1 deg.
particles = danton.particles(10, energy=1E+09, elevation=1.0)

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Inspect the results.
print(result)
