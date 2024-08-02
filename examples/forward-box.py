#! /usr/bin/env python3
import danton

# Create a Danton simulation interface.
simulation = danton.Simulation()

# Configure the Monte Carlo for the forward sampling of tau decays (using the
# longitudinal approximlation).
simulation.mode = "forward"
simulation.tau_decays = True
simulation.longitudinal = True

# Set a spherical Earth (using the PREM model).
simulation.geometry.geodesic = "PREM"

# Define a local box-area for collecting tau decays.
box = simulation.box(size=[1E+05, 1E+05, 1E+04])

# Generate a bunch of neutrinos targeting the box (and preselect -inclusively-
# events that might lead to decays with an elevation angle within [-1, 1] deg).
particles, size = simulation.particles() \
    .powerlaw(1E+07, 1e+11) \
    .solid_angle(elevation=[-1.0, 1.0]) \
    .target(box) \
    .generate(1000)

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Inspect the results.
secondaries = result.secondaries
inside = box.inside(secondaries) # Select decays occuring inside the box.
print(box.to_local(secondaries)) # Print secondaries coordinates in box frame.
