"""
Created on Thu May 25 10:21:52 2017
@author: Rachael A. Mansbach
A simple HOOMD script to test rigid body particles with patches on them, mainly for timing purposes
Three overlapping beads with excluded volume interactions plus one LJ particle on each end of the same size plus six patches above and below the body of the particle       

Put 1000 molecules into a box, 500 of type A and 500 of type B
For the time being, the only difference between type A and type B is their SC
radius and SC well depth

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

def make_DXXX_monomer(bbbtype,bbbtype2,scbtype,abtype,beadR,ebeadR,lbeadR,sticky_theta,
                      rigid):
    """
    Updates rigid body to define a particular bead type in the correct shape
    
    ----------
    Parameters
    ----------
    bbbtype: string
        BB bead type name
    bbbtype2: string
        BB bead 2 type name
    scbtype: string
        SC bead type name
    abtype: string
        A bead type name
    beadR: float
        bead radius
    ebeadR: float
        offset bead radius
    lbeadR: float
        little bead radius
    sticky_theta: float
        angle to place little beads at
    rigid: HOOMD rigid object
    
    -------
    Returns
    ------
    Updated rigid
    """
    rigid.set_param(bbbtype,positions=[(-beadR,0,0),(beadR,0,0),(-2*beadR-ebeadR,0,0),
                               (2*beadR+ebeadR,0,0),
                               (0,(beadR-lbeadR)*np.cos(sticky_theta),
                                (beadR-lbeadR)*np.sin(sticky_theta)),
                                (0,-(beadR-lbeadR)*np.cos(sticky_theta),
                                 (beadR-lbeadR)*np.sin(sticky_theta)),
                                 (beadR,(beadR-lbeadR)*np.cos(sticky_theta),
                                  (beadR-lbeadR)*np.sin(sticky_theta)),
                                  (beadR,-(beadR-lbeadR)*np.cos(sticky_theta),
                                   (beadR-lbeadR)*np.sin(sticky_theta)),
                                   (-beadR,(beadR-lbeadR)*np.cos(sticky_theta),
                                    (beadR-lbeadR)*np.sin(sticky_theta)),
                                    (-beadR,
                                     -(beadR-lbeadR)*np.cos(sticky_theta),
                                     (beadR-lbeadR)*np.sin(sticky_theta)),
                                     (0,(beadR-lbeadR)*np.cos(sticky_theta),
                                      -(beadR-lbeadR)*np.sin(sticky_theta)),
                                      (0,-(beadR-lbeadR)*np.cos(sticky_theta),
                                       -(beadR-lbeadR)*np.sin(sticky_theta)),
                                       (beadR,
                                        (beadR-lbeadR)*np.cos(sticky_theta),
                                        -(beadR-lbeadR)*np.sin(sticky_theta)),
                                        (beadR,
                                         -(beadR-lbeadR)*np.cos(sticky_theta),
                                         -(beadR-lbeadR)*np.sin(sticky_theta)),
                                         (-beadR,
                                          (beadR-lbeadR)*np.cos(sticky_theta),
                                          -(beadR-lbeadR)*np.sin(sticky_theta)),
                                          (-beadR,
                                           -(beadR-lbeadR)*np.cos(sticky_theta),
                                           -(beadR-lbeadR)*np.sin(sticky_theta))],
                                           types=[bbbtype2,bbbtype2,scbtype,scbtype,abtype,
                                                  abtype,abtype,abtype,abtype,abtype,
                                                  abtype,abtype,abtype,abtype,abtype,
                                                  abtype]) #type of particle has to be different from central particle because otherwise they are defined recursively
    
    return rigid
        
beadMass = 1
beadR = 0.5
ebeadR=[1.5/(2.**(5./6.)),1.0/(2.**(5./6.))]
lbeadR = 0.125
sticky_theta = 10.*(np.pi/180.) #angle down the main sphere of the little spheres
seed = 123452
sigslj = 2*beadR #shifted LJ sigma is basically particle diameter
rcutslj = sigslj
rcutlb = 2.0
siglb =[1.5,1.0] #diameter of particle is roughly the LJ min, rm = 2^(1/6) * sigma
sigls = 2*lbeadR / 2.**(1./6.)
siglbb = 2*beadR / 2.**(1./6.)
#initial rough estimates from Martini oriented dimer PMF's, lseps = 10.0, lbeps = 4.5
lseps=2.5
lbeps=[0.2,0.9]
lbb=1.0
offset=0.025

context.initialize()
Ntypes = 2
latr = 19.4945
a1 = [latr,0,0]
a2 = [0,latr,0]
a3 = [0,0,latr]
pos =[[0.1*latr,0.33*latr,0.33*latr],[0.1*latr,0.66*latr,0.33*latr],
      [0.1*latr,0.33*latr,0.75*latr],[0.1*latr,0.66*latr,0.66*latr],
      [0.9*latr,0.33*latr,0.5*latr],[0.9*latr,0.66*latr,0.5*latr],
      [0.9*latr,0.33*latr,0.66*latr],[0.9*latr,0.66*latr,0.75*latr],
      [0.33*latr,0.1*latr,0.33*latr],[0.66*latr,0.1*latr,0.33*latr],
      [0.33*latr,0.1*latr,0.75*latr],[0.66*latr,0.1*latr,0.75*latr],
      [0.33*latr,0.9*latr,0.5*latr],[0.66*latr,0.9*latr,0.5*latr],
      [0.33*latr,0.9*latr,0.66*latr],[0.66*latr,0.9*latr,0.66*latr],
      [0.25*latr,0.25*latr,0.25*latr],[0.25*latr,0.25*latr,0.75*latr],
      [0.75*latr,0.25*latr,0.25*latr],[0.75*latr,0.25*latr,0.75*latr]]
types = ['A','B']
mtypes = ["EB","EB","EB","EB","EB","EB","EB","EB","EA","EB","EB","EB","EB","EA","EB","EB","EA","EA","EA","EB"]
ucell = hoomd.lattice.unitcell(20, a1, a2, a3, dimensions=3, position=pos,
                               type_name=mtypes)

system = init.create_lattice(unitcell = ucell, n=[8,8,8]) 
rigid = md.constrain.rigid()
ind = 0
system.particles.types.add('E2')
system.particles.types.add('LS')
for t in types:
    system.particles.types.add('LB'+t)
    lbeadR-=offset
    rigid = make_DXXX_monomer('E'+t,'E2','LB'+t,'LS',beadR,ebeadR[ind],lbeadR,
                              sticky_theta,rigid)
    lbeadR+=offset
    ind += 1
    
rigid.create_bodies()

groupR = group.rigid_center()
groupLBs = [group.type('LB'+t) for t in types]
groupLS = group.type('LS')
groupBB = (group.type('EA') or group.type('E2') or group.type('EB'))

for p in groupR:
    #set moment of inertia based on overlapping spheres and radius = 0.5 with rigid beads laid out along x axis
    Ixx = 3*isphere(beadMass,beadR,0)+2*isphere(3.75*beadMass,beadR,0)
    Iyy = isphere(beadMass,beadR,0) + 2*isphere(beadMass,beadR,beadR) + 2*isphere(3.75*beadMass,beadR,3*beadR)
    Izz = Iyy
    p.moment_inertia = (Ixx,Iyy,Izz)
 #set masses and diameters of beads
for p in groupBB: 
    p.mass = beadMass
    p.diameter = 2*beadR
for i in range(len(groupLBs)):    
    for p in groupLBs[i]:
        p.mass = 3.75*beadMass
        p.diameter = 2*ebeadR[i]

for p in groupLS:
    p.mass = 0.0
    p.diameter = 2*lbeadR
#initialize md neighbor list and simple force field with slightly soft repulsion for repulsive beads (remember you can't have discontinuous potentials in MD!)
#eventually this should be done with a function and some parameter input (see Andy L's code for inspiration)
epsDictSLJ = dict()
for ti in range(len(types)):
    t1 = types[ti]
    epsDictSLJ[('E'+t1,'E'+t1)] = 1.0
    epsDictSLJ[('E'+t1,'E2')] = 1.0
    epsDictSLJ[('E2','LB'+t)] = 1.0
    for tj in range(ti,len(types)):
        t2 = types[tj]
        epsDictSLJ[('LB'+t1,'LB'+t2)] = 1.0
    for tk in range(len(types)):
        t3 = types[tk]      
        epsDictSLJ[('E'+t1,'LB'+t3)] = 1.0
epsDictSLJ[('E2','E2')] = 1.0



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
types = ['EA','EB','E2','LS','LBA','LBB']
#don't need to turn on slj if lj is taking care of it
if (lbb==0) or (lbeps[0]==0) or (lbeps[1]==0):
	slj = rd.setupSLJ(types,epsDictSLJ,sigDictSLJ,rcutslj,nl,2*beadR)


epsDictLJ = dict()
epsDictLJ[('LS','LS')] = lseps
epsDictLJ[('LBA','LBA')] = lbeps[0]
epsDictLJ[('LBB','LBB')] = lbeps[1]
epsDictLJ[('LBA','LBB')] = MIXE(lbeps[0],lbeps[1])
epsDictLJ[('LS','LBA')] = MIXE(lseps,lbeps[0])
epsDictLJ[('LS','LBB')] = MIXE(lseps,lbeps[1])
epsDictLJ[('EA','EA')] = lbb
epsDictLJ[('EA','E2')] = lbb
epsDictLJ[('E2','E2')] = lbb
epsDictLJ[('EA','EB')] = lbb
epsDictLJ[('EB','EB')] = lbb
epsDictLJ[('EB','E2')] = lbb
epsDictLJ[('EA','LS')] = MIXE(lbb,lseps)
epsDictLJ[('EB','LS')] = MIXE(lbb,lseps)
epsDictLJ[('E2','LS')] = MIXE(lbb,lseps)
epsDictLJ[('EA','LBA')] = MIXE(lbb,lbeps[0])
epsDictLJ[('EA','LBB')] = MIXE(lbb,lbeps[1])
epsDictLJ[('EB','LBA')] = MIXE(lbb,lbeps[0])
epsDictLJ[('EB','LBB')] = MIXE(lbb,lbeps[1])
epsDictLJ[('E2','LBA')] = MIXE(lbb,lbeps[0])
epsDictLJ[('E2','LBB')] = MIXE(lbb,lbeps[1])

sigDictLJ = dict()
sigDictLJ[('LS','LS')] = sigls
sigDictLJ[('LBA','LBA')] = siglb[0]
sigDictLJ[('LBB','LBB')] = siglb[1]
sigDictLJ[('LS','LBA')] = MIXS(sigls,siglb[0])
sigDictLJ[('LS','LBB')] = MIXS(sigls,siglb[1])
sigDictLJ[('LBA','LBB')] = MIXS(siglb[0],siglb[1])
sigDictLJ[('EA','EA')] = siglbb
sigDictLJ[('EB','EB')] = siglbb
sigDictLJ[('EA','EB')] = siglbb
sigDictLJ[('E2','E2')] = siglbb
sigDictLJ[('EA','E2')] = siglbb
sigDictLJ[('EB','E2')] = siglbb
sigDictLJ[('EA','LBA')] = MIXS(siglbb,siglb[0])
sigDictLJ[('EA','LBB')] = MIXS(siglbb,siglb[1])
sigDictLJ[('EB','LBA')] = MIXS(siglbb,siglb[0])
sigDictLJ[('EB','LBB')] = MIXS(siglbb,siglb[1])
sigDictLJ[('E2','LBA')] = MIXS(siglbb,siglb[0])
sigDictLJ[('E2','LBB')] = MIXS(siglbb,siglb[1])
sigDictLJ[('EA','LS')] = MIXS(siglbb,sigls)
sigDictLJ[('EB','LS')] = MIXS(siglbb,sigls)
sigDictLJ[('E2','LS')] = MIXS(siglbb,sigls)

lj = rd.setupLJ(types,epsDictLJ,sigDictLJ,rcutlb,nl)

gsd_restart = dump.gsd(filename="restart_run2.gsd",group=group.all(),truncate=True,period=200000,phase=0)
dump.gsd(filename='mols10000_250-02-09-150-100_run2.gsd',period=5000,group=group.all(),phase=0)
dump.gsd(filename='mols10000_250-02-09-150-100_short_run2.gsd',period=25000,group=group.all(),phase=0)
md.integrate.mode_standard(dt=1e-3, aniso=True)
langevin = md.integrate.langevin(group=groupR, kT=1.0, seed=seed,dscale=1.0,tally=True)
analyze.log(filename='mols10000_250-02-09-150-100_run2.log',quantities=['potential_energy','kinetic_energy','langevin_reservoir_energy_rigid_center','temperature','pressure'],period=50000,phase=0)
try:
	run_upto(2e2,limit_multiple=200000)
except WalltimeLimitReached:
	pass
gsd_restart.write_restart()
if (lbb==0) or (lbeps[0]==0) or (lbeps[1]==0):
    slj.disable()
lj.disable()
langevin.disable()
