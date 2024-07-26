#! /usr/bin/env python3
import danton

simulation = danton.Simulation()
simulation.mode = "forward"
simulation.geometry.geodesic = "PREM"
particles = danton.particles(100, pid=16, energy=1E+10, elevation=1.0)
result = simulation.run(particles)
print(result)
