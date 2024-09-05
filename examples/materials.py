#! /usr/bin/env python3
import danton

# Instanciate a simulation with custom materials.
simulation = danton.Simulation(
    materials = "examples/rocks.toml", # See the corresponding file.
)

# By default, the topography layer is made of (standard) "Rock".
print(f"topography material: {simulation.geometry.material}")

# Let us change this to Limestone (defined in examples/rocks.toml).
simulation.geometry.material = "Limestone"

# The topography bulk density is `None` by default, meaning that the native
# material density is used.
print(f"topography density: {simulation.geometry.density}")

# Let us set a lower bulk density, e.g. to account for porous rocks.
simulation.geometry.density = 2E+03 # kg / m3
