#! /usr/bin/env python3
import danton


# Simulation settings.
energy = 1E+10 # eV
latitude, longitude = -30.705, -68.263 # deg
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

# Convert the ECEF direction of incoming neutrinos to horizontal coordinates at
# the box centre.
ecef_position = simulation.geometry.to_ecef(
    latitude = latitude,
    longitude = longitude
)
coordinates = simulation.geometry.from_ecef(ecef_position, ecef_direction)
azimuth, elevation = coordinates["azimuth"], coordinates["elevation"]

print(f"azimuth = {azimuth:.2f} deg, elevation = {elevation:.2f} deg")

# Generate a bunch of neutrinos targeting the box.
particles = simulation.particles() \
    .energy(energy) \
    .direction(azimuth=azimuth, elevation=elevation) \
    .target(box) \
    .generate(1000)

area = box.projected_surface(azimuth=azimuth, elevation=elevation)
print(f"generation area = {area * 1E-06:.3f} km^2")

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Inspect the results.
secondaries = result.secondaries
inside = box.inside(secondaries) # Select decays occuring inside the box.
print(simulation.geometry.to_ecef(secondaries)) # Print secondaries coordinates
                                                # in ECEF.
