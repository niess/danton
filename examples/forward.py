#! /usr/bin/env python3
import danton

# Create a Danton simulation interface.
simulation = danton.Simulation()

# Configure the Monte Carlo for forward flux sampling (using the longitudinal
# approximlation).
simulation.mode = "forward"
simulation.tau_decays = False
simulation.longitudinal = True

# Set a spherical Earth (using PREM model).
simulation.geometry.geodesic = "PREM"

# Generate a bunch of neutrinos (that would emerge with an elevation angle of 1
# deg at the equator).
particles = danton.particles(100, pid=16, energy=1E+10, elevation=1.0)

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Inspect the results.
print(result)
