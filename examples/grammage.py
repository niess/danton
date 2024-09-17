#! /usr/bin/env python3
import danton

simulation = danton.Simulation()
simulation.mode = "grammage"
simulation.record_steps = True
simulation.geometry.geoid = "PREM81"
particles = danton.particles(elevation=-30.0)
result = simulation.run(particles)
print(result.steps)
