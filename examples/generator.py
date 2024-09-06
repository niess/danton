#! /usr/bin/env python3
import danton

N = 1000
simulation = danton.Simulation()
particles = simulation.particles()       \
    .powerlaw(1E+06, 1E+12, exponent=-2) \
    .solid_angle(elevation=[-1, 1])      \
    .generate(N)

print(particles)
