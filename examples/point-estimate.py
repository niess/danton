#! /usr/bin/env python3
import danton

def flux(energy):
    """Some neutrino flux."""
    return 1.5E-04 # 1 / (GeV m^2 sr s)

N = 1000
simulation = danton.Simulation(mode="backward")
taus = danton.particles(N, # Create N identical particles.
    energy = 1E+10, # GeV
    altitude = 100, # m
    elevation = 1   # deg
)
primaries = simulation.run(taus).primaries # Backward propagate the taus.
density = sum(primaries["weight"] * flux(primaries["energy"])) / N # Weighted average.

print(f"density = {density:.3E} 1 / (GeV m^3 sr s)")
