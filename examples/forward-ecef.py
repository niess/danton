#! /usr/bin/env python3
import danton


# Simulation settings.
energy = 1E+10 # eV
latitude, longitude = -30.705, -68.263 # deg
altitude = 700.0 # m
ecef_direction = [0.96339066, 0.19163032, 0.1875]

# Create a Danton simulation interface.
simulation = danton.Simulation()

# Configure the Monte Carlo for the forward sampling of tau decays (using the
# longitudinal approximlation).
simulation.mode = "forward"
simulation.tau_decays = True
simulation.longitudinal = True

# Set an elliptical Earth (using the WGS84 ellipsoid).
simulation.geometry.geoid = "WGS84"

# Define a local box-area for collecting tau decays.
box = simulation.box(
    latitude = latitude,
    longitude = longitude,
    size = [1E+05, 1E+05, 1E+04]
)

# Generate a bunch of neutrinos targeting the box.
particles = simulation.particles() \
    .energy(energy) \
    .direction(ecef_direction) \
    .target(box) \
    .generate(1000)

area = box.projected_area(direction=ecef_direction)
print(f"generation area = {area * 1E-06:.3f} km^2")

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Inspect the results.
secondaries = result.secondaries
inside = box.inside(secondaries) # Select decays occuring inside the box.
print(simulation.geometry.to_ecef(secondaries)) # Print secondaries coordinates
                                                # in ECEF.
