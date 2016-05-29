## lj.py: a Monte Carlo simulation of the Lennard-Jones fluid, in the NVT ensemble.
import random
import sys
import math

# Simulation Parameters
N       = 500
sigma   = 1
epsilon = 1
trunc   = 3*sigma
truncsq = trunc**2
equil   = int(sys.argv[1])
prod    = 2.5e8
temp    = 8.5e-1
density = 1.0e-3
L =  (N/density)**(1.0/3.0)
halfL = L/2

particles = []

# Some helper functions
def wrap(particle):
	this_particle = particle
	if particle[0] > L:
		this_particle[0] -= L
	elif particle[0] < 0:
		this_particle[0] += L
	if particle[1] > L:
		this_particle[1] -= L
	elif particle[1] < 0:
		this_particle[1] += L
	if particle[2] > L:
		this_particle[2] -= L
	elif particle[2] < 0:
		this_particle[2] += L	
	return this_particle

def distancesq(particle1, particle2):
	'''Gets the squared distance between two particles, applying the minimum image convention.'''
	# Calculate distances
	dx = particle1[0]-particle2[0]
	dy = particle1[1]-particle2[1]
	dz = particle1[2]-particle2[2]

	# Minimum image convention
	if dx > halfL:
		dx -= L
	elif dx < -halfL:
		dx += L
	if dy > halfL:
		dy -= L
	elif dy < -halfL:
		dy += L
	if dz > halfL:
		dz -= L
	elif dz < -halfL:
		dz += L

	return dx**2+dy**2+dz**2

def energy(particles):
	'''Gets the energy of the system'''
	energy = 0
	for particle1 in range(0, len(particles)-1):
		for particle2 in range(particle1+1, len(particles)):
			dist = distancesq(particles[particle1], particles[particle2])
			if dist <= truncsq:
				energy += (1/dist**6)-(1/dist**3)
	return energy

def particleEnergy(particle, particles, p):
	part_energy = 0
	i = 0
	for particle2 in particles:
		if i != p:
			dist = distancesq(particle, particle2)
			if dist <= truncsq:
				part_energy += (1/dist**6)-(1/dist**3)
		i += 1
	return part_energy

def writeEnergy(step, en):
	with open('energy', 'a') as f:
		f.write('{0} {1}\n'.format(step, en))

# Pack the box:
for particle in range(0, N):
	x_coord = random.uniform(0, L)
	y_coord = random.uniform(0, L)
	z_coord = random.uniform(0, L)
	particles.append([x_coord, y_coord, z_coord])

# Calculate initial energy
en = 0
en = energy(particles)

print('Initial energy: {0}'.format(en))

# MC
for step in range(0, equil):
	sys.stdout.write("\rSTEP: {0} Energy: {1}".format(step, en))
	sys.stdout.flush()
	p = 0
	for particle in particles:
		this_particle = particle
		prev_E = particleEnergy(particle, particles, p)
		this_particle[0] += random.uniform(-1, 1)
		this_particle[1] += random.uniform(-1, 1)
		this_particle[2] += random.uniform(-1, 1)
		this_particle = wrap(this_particle)
		new_E = particleEnergy(this_particle, particles, p)
		deltaE = new_E - prev_E
		# Check energy
		if deltaE < 0:
			particles[p] = this_particle
			en += new_E - prev_E
		else:
			rand = random.random()
			if math.exp(-deltaE/temp) > rand:
				particles[p] = this_particle
				en += new_E-prev_E
		p += 1
	writeEnergy(str(step), str(en))
