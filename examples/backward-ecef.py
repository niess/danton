#! /usr/bin/env python3
import danton


# Simulation settings.
energy = [1E+06, 1E+12] # GeV
latitude, longitude = -30.705, -68.263 # deg
altitude = 700.0 # m
ecef_direction = [0.96339066, 0.19163032, 0.1875]

# Create a Danton simulation interface.
simulation = danton.Simulation()

# Configure the Monte Carlo for the backward sampling of tau decays (using the
# longitudinal approximlation).
simulation.mode = "backward"
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

# Generate a bunch of tau decays inside the box.
particles, N = simulation.particles()  \
    .powerlaw(*energy, exponent=-1)    \
    .direction(ecef_direction)         \
    .inside(box)                       \
    .generate(10000)

# Run the Monte Carlo simulation.
result = simulation.run(particles)

# Apply a primary flux model.
def flux(energy):
    """Some neutrino flux."""
    return 1.5E-04 / energy**2 # 1 / (GeV m^2 sr s)

primaries = result.primaries
primaries = primaries[primaries["pid"] == 16] # Select tau neutrinos.
primaries["weight"] *= flux(primaries["energy"])

# Compute (a MC estimate of) the rate of tau decays over the box volume.
rate = sum(primaries["weight"]) / N
sigma = ((sum(primaries["weight"]**2) / N - rate**2) / N)**0.5
print(f"rate = {rate:.3E} +- {sigma:.3E}  [Hz / sr]")
