"""
Created on Thu May 25 10:21:52 2017
@author: Rachael A. Mansbach
A simple HOOMD script to test rigid body particles with patches on them, mainly for timing purposes
Three overlapping beads with excluded volume interactions plus one LJ particle on each end of the same size plus six patches above and below the body of the particle       

Put 378 molecules into a box of size 24.3 x 24.3 x 24.3 nm (after Martini simulations at highest concentration)

Units: 
distance d* = 1 nm
mass m* = 102 amu from the mass of one of the aromatic ring residues (8 carbons + 6 hydrogens)
energy eps* = kb T at 298 K (16 kb T ~ dimer FE change)
temp T = 1
resulting t* = sqrt((m* d*^2)/eps*) = 6.4 ps
"""
import hoomd,imp
from hoomd import *
from hoomd import md
import numpy as np

import gsd.hoomd
rd = imp.load_source('rd','rundebughelpers.py')
#function to determine the moment of inertia of a sphere at a distance of r from the origin, with a mass of m and a radius R
def isphere(m,R,r):
    I = (2./5.)*m*R*R+m*r*r
    return I
    
#function to determine the moment of inertia of a line of spheres spaced symmetrically about the origin, with masses m and radii R, and a spacing r
#n is the number of spheres to each side of the central sphere, so the total number is 2n+1
def isphereline(m,R,r,n):
    I = isphere(m,R,0)
    for i in range(n):
        I += 2*isphere(m,R,(i+1)*r)
    return I    
    
# Lorentz-Bertholot Mixing Rules
def MIXE(eps_i, eps_j):
	return np.sqrt(eps_i*eps_j)

def MIXS(sigma_i, sigma_j):
	return (sigma_i+sigma_j)/2.0

        
beadMass = 1
beadR = 0.5
ebeadR=ERAD/(2.**(5./6.))
lbeadR = 0.125
sticky_theta = 10.*(np.pi/180.) #angle down the main sphere of the little spheres
seed = SSSRUN
sigslj = 2*beadR #shifted LJ sigma is basically particle diameter
rcutslj = sigslj
rcutlb = 2.0
siglb =ERAD #diameter of particle is roughly the LJ min, rm = 2^(1/6) * sigma
sigls = 2*lbeadR / 2.**(1./6.)
siglbb = 2*beadR / 2.**(1./6.)
#initial rough estimates from Martini oriented dimer PMF's, lseps = 10.0, lbeps = 4.5
lseps=AAA
lbeps=SCSCSC
lbb=1.0
offset=0.025

context.initialize()
system = init.create_lattice(unitcell = hoomd.lattice.sc(a=7.1814, type_name='E'), n=[22,22,22]) #replicate a lattice of types 7 x 6 x 9 times for a total of 378 molecules
#(rigid,groupR,groupLB,groupLS) = rd.rigid5init(system,beadMass,beadR,lbeadR,sticky_theta)
#system = init.read_gsd(filename='/scratch1/scratchdirs/mansbac2/patchy/rigid_assembly/offset_little/init.gsd',time_step=0)
#(rigid,groupR,groupLB,groupLS) = rd.rigid5init(system,beadMass,beadR,lbeadR,sticky_theta)
#add secondary particle types to simulation
system.particles.types.add('E2')
system.particles.types.add('LB')
system.particles.types.add('LS')
rigid = md.constrain.rigid()
lbeadR-=offset
rigid.set_param('E',positions=[(-beadR,0,0),(beadR,0,0),(-2*beadR-ebeadR,0,0),(2*beadR+ebeadR,0,0),(0,(beadR-lbeadR)*np.cos(sticky_theta),(beadR-lbeadR)*np.sin(sticky_theta)),(0,-(beadR-lbeadR)*np.cos(sticky_theta),(beadR-lbeadR)*np.sin(sticky_theta)),(beadR,(beadR-lbeadR)*np.cos(sticky_theta),(beadR-lbeadR)*np.sin(sticky_theta)),(beadR,-(beadR-lbeadR)*np.cos(sticky_theta),(beadR-lbeadR)*np.sin(sticky_theta)),(-beadR,(beadR-lbeadR)*np.cos(sticky_theta),(beadR-lbeadR)*np.sin(sticky_theta)),(-beadR,-(beadR-lbeadR)*np.cos(sticky_theta),(beadR-lbeadR)*np.sin(sticky_theta)),(0,(beadR-lbeadR)*np.cos(sticky_theta),-(beadR-lbeadR)*np.sin(sticky_theta)),(0,-(beadR-lbeadR)*np.cos(sticky_theta),-(beadR-lbeadR)*np.sin(sticky_theta)),(beadR,(beadR-lbeadR)*np.cos(sticky_theta),-(beadR-lbeadR)*np.sin(sticky_theta)),(beadR,-(beadR-lbeadR)*np.cos(sticky_theta),-(beadR-lbeadR)*np.sin(sticky_theta)),(-beadR,(beadR-lbeadR)*np.cos(sticky_theta),-(beadR-lbeadR)*np.sin(sticky_theta)),(-beadR,-(beadR-lbeadR)*np.cos(sticky_theta),-(beadR-lbeadR)*np.sin(sticky_theta))],types=['E2','E2','LB','LB','LS','LS','LS','LS','LS','LS','LS','LS','LS','LS','LS','LS']) #type of particle has to be different from central particle because otherwise they are defined recursively
lbeadR+=offset
rigid.create_bodies()

groupR = group.rigid_center()
groupLB = group.type('LB')
groupLS = group.type('LS')
groupBB = (group.type('E') or group.type('E2'))

for p in group.type('E'):
    #set moment of inertia based on overlapping spheres and radius = 0.5 with rigid beads laid out along x axis
    Ixx = 3*isphere(beadMass,beadR,0)+2*isphere(3.75*beadMass,beadR,0)
    Iyy = isphere(beadMass,beadR,0) + 2*isphere(beadMass,beadR,beadR) + 2*isphere(3.75*beadMass,beadR,3*beadR)
    Izz = Iyy
    p.moment_inertia = (Ixx,Iyy,Izz)
 #set masses and diameters of beads
for p in groupBB: 
    p.mass = beadMass
    p.diameter = 2*beadR
    
for p in groupLB:
    p.mass = 3.75*beadMass
    p.diameter = 2*ebeadR

for p in groupLS:
    p.mass = 0.0
    p.diameter = 2*lbeadR
#initialize md neighbor list and simple force field with slightly soft repulsion for repulsive beads (remember you can't have discontinuous potentials in MD!)
#eventually this should be done with a function and some parameter input (see Andy L's code for inspiration)
epsDictSLJ = dict()
epsDictSLJ[('E','E')] = 1.0
epsDictSLJ[('E2','E2')] = 1.0
epsDictSLJ[('E','E2')] = 1.0
epsDictSLJ[('E','LB')] = 1.0
epsDictSLJ[('E2','LB')] = 1.0
epsDictSLJ[('LB','LB')] = 1.0

sigDictSLJ = dict()
if lbb == 0:
	sigDictSLJ[('E','E')] = sigslj
	sigDictSLJ[('E2','E2')] = sigslj
	sigDictSLJ[('E','E2')] = sigslj
	sigDictSLJ[('E','LB')] = sigslj
	sigDictSLJ[('E2','LB')] = sigslj
else:
	sigDictSLJ[('E','E')] = sigslj
	sigDictSLJ[('E2','E2')] = sigslj
	sigDictSLJ[('E','E2')] = sigslj
	sigDictSLJ[('E','LB')] = sigslj
	sigDictSLJ[('E2','LB')] = sigslj

if lbeps == 0:
	sigDictSLJ[('LB','LB')] = sigslj
else:
	sigDictSLJ[('LB','LB')] = 0.0

nl = md.nlist.tree(r_buff=0.6)
types = ['E','E2','LS','LB']
#don't need to turn on slj if lj is taking care of it
if (lbb==0) or (lbeps==0):
	slj = rd.setupSLJ(types,epsDictSLJ,sigDictSLJ,rcutslj,nl,2*beadR)


epsDictLJ = dict()
epsDictLJ[('LS','LS')] = lseps
epsDictLJ[('LB','LB')] = lbeps
epsDictLJ[('LS','LB')] = MIXE(lseps,lbeps)
epsDictLJ[('E','E')] = lbb
epsDictLJ[('E','E2')] = lbb
epsDictLJ[('E2','E2')] = lbb
epsDictLJ[('E','LS')] = MIXE(lbb,lseps)
epsDictLJ[('E2','LS')] = MIXE(lbb,lseps)
epsDictLJ[('E','LB')] = MIXE(lbb,lbeps)
epsDictLJ[('E2','LB')] = MIXE(lbb,lbeps)

sigDictLJ = dict()
sigDictLJ[('LS','LS')] = sigls
sigDictLJ[('LB','LB')] = siglb
sigDictLJ[('LS','LB')] = MIXS(sigls,siglb)
sigDictLJ[('E','E')] = siglbb
sigDictLJ[('E2','E2')] = siglbb
sigDictLJ[('E','E2')] = siglbb
sigDictLJ[('E','LB')] = MIXS(siglbb,siglb)
sigDictLJ[('E2','LB')] = MIXS(siglbb,siglb)
sigDictLJ[('E','LS')] = MIXS(siglbb,sigls)
sigDictLJ[('E2','LS')] = MIXS(siglbb,sigls)

lj = rd.setupLJ(types,epsDictLJ,sigDictLJ,rcutlb,nl)
gsd_restart = dump.gsd(filename="restart_runRUN.gsd",group=group.all(),truncate=True,period=200000,phase=0)
dump.gsd(filename='mols10000_AAA-SCSCSC-ERAD_runRUN.gsd',period=5000,group=group.all(),phase=0)
dump.gsd(filename='mols10000_AAA-SCSCSC-ERAD_short_runRUN.gsd',period=25000,group=group.all(),phase=0)
md.integrate.mode_standard(dt=1e-3, aniso=True)
langevin = md.integrate.langevin(group=groupR, kT=1.0, seed=seed,dscale=1.0,tally=True)
analyze.log(filename='mols10000_AAA-SCSCSC-ERAD_runRUN.log',quantities=['potential_energy','kinetic_energy','langevin_reservoir_energy_rigid_center','temperature','pressure'],period=50000,phase=0)
try:
	run_upto(2e7,limit_multiple=200000)
except WalltimeLimitReached:
	pass
gsd_restart.write_restart()
if (lbb==0) or (lbeps==0):
    slj.disable()
lj.disable()
langevin.disable()
