#! /usr/bin/env python3
import danton

# Create a Danton simulation interface.
simulation = danton.Simulation()

# Configure the Monte Carlo for the backward sampling of tau decays 
simulation.mode = "backward"
simulation.tau_decays = True

# Set a spherical Earth (using the PREM model).
simulation.geometry.geoid = "PREM81"

# Define a local box-area for collecting tau decays.
box = simulation.box(size=[1E+05, 1E+05, 1E+04])

# Generate a bunch of tau decays within the box volume (applying a preselction).
particles, size = simulation.particles() \
    .powerlaw(1E+06, 1E+12) \
    .solid_angle(elevation=[-1.0, 1.0]) \
    .inside(box, limit=True) \
    .generate(100)

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Inspect the results.
primaries = result.primaries
print(primaries)
