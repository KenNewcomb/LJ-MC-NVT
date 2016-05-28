## lj.py: a Monte Carlo simulation of the Lennard-Jones fluid, in the NVT ensemble.
import random

# Simulation Parameters
N       = 500
sigma   = 1
epsilon = 1
trunc   = 3*sigma 
equil   = 5e7
prod    = 2.5e8
temp    = 8.50e-1
density = 1.0e-3

x_coords = []
y_coords = []
z_coords = []

# Pack the box:
for particle in range(0, N):
	x_coords.append(random.uniform(0, (N/density)**(1.0/3.0)))
	y_coords.append(random.uniform(0, (N/density)**(1.0/3.0)))
	z_coords.append(random.uniform(0, (N/density)**(1.0/3.0)))

# MC
