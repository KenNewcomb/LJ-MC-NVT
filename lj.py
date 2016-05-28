## lj.py: a Monte Carlo simulation of the Lennard-Jones fluid, in the NVT ensemble.
import random

# Simulation Parameters
N       = 500
sigma   = 1
epsilon = 1
trunc   = 3*sigma 
equil   = int(5e7)
prod    = 2.5e8
temp    = 8.50e-1
density = sigma**3

particles = []

# Some helper functions
def wrap(particle):
	if particle[0] > sigma:
		particle[0] -= sigma
	elif particle[0] < 0:
		particle[0] += sigma
	if particle[1] > sigma:
		particle[1] -= sigma
	elif particle[1] < 0:
		particle[1] += sigma
	if particle[2] > sigma:
		particle[2] -= sigma
	elif particle[2] < 0:
		particle[2] += sigma
	
	return particle
def distance(particle1, particle2):
	'''Gets the distance between two particles, applying the minimum image convention.'''
	# Calculate distances
	dx = particle1[0]-particle2[0]
	dy = particle1[1]-particle2[1]
	dz = particle1[2]-particle2[2]

	# Minimum image convention
	if dx > sigma/2:
		dx -= sigma
	elif dx < -sigma/2:
		dx += sigma
	if dy > sigma/2:
		dy -= sigma
	elif dy < -sigma/2:
		dy += sigma
	if dz > sigma/2:
		dz -= sigma
	elif dz < -sigma/2:
		dz += sigma

	return (dx**2+dy**2+dz**2)**(0.5)

def energy(particles):
	'''Gets the energy of the system'''
	energy = 0
	for particle1 in range(0, len(particles)-1):
		for particle2 in range(particle1+1, len(particles)):
			dist = distance(particles[particle1], particles[particle2])
			if dist <= trunc:
				energy += epsilon*((1/dist**12)-(1/dist**6))
	return energy

# Pack the box:
for particle in range(0, N):
	x_coord = random.uniform(0, (N/density)**(1.0/3.0))
	y_coord = random.uniform(0, (N/density)**(1.0/3.0))
	z_coord = random.uniform(0, (N/density)**(1.0/3.0))
	particles.append([x_coord, y_coord, z_coord])

# Calculate initial energy
en = 0
en = energy(particles)

print('Initial energy: {0}'.format(en))

# MC
for step in range(0, equil):
	print("STEP: {0}".format(step))
	p = 0
	for particle in particles:
		this_particle = particle
		this_particle[0] += random.uniform(-0.1, 0.1)
		this_particle[1] += random.uniform(-0.1, 0.1)
		this_particle[2] += random.uniform(-0.1, 0.1)
		
		this_particle = wrap(this_particle)

		particles[p] = this_particle

		# Check energy
		if energy(particles) < en:
			en = energy(particles)
			print("NEW ENERGY: {0}".format(en))
		else:
			particles[p] = particle
		p += 1
