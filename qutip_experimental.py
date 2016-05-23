from qutip import *
from pylab import *
from scipy import constants as cp
from numpy import random


def qubit_state(theta,phi):
    return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])


def circ_motion(t, args):
    r=args['cOffset'] + array( [ -args['D']/2. + args['D']/sqrt(2)*sin(2*pi/args['tau'] * t)  ,  -args['D']/2. + args['D']/sqrt(2)*cos(2*pi/args['tau'] * t)  ,  args['d'] ] )
    #allow random path jitter (simulates not perfect movement of mems stage) ~~1nm?
    if args['pJit']:
        # add gavins function here!
        r+=random.normal(0, args['rstd'], 3)
        # use below for test plot
        #r=args['cOffset'] + array( [ -args['D']/2. + random.normal(args['D']/sqrt(2)*sin(2*pi/args['tau'] * t), args['rstd'])  ,  -args['D']/2. + random.normal(args['D']/sqrt(2)*cos(2*pi/args['tau'] * t), args['rstd'])  ,  random.normal(args['d'], args['rstd']) ] )
    return r

def H_RWA(t, args):
    r=args['r']
    phase=exp(4j*args['Delta'])
    if args['circ']:
        #implements circular motion from with radius D/sqrt2 where D is the separation between donors. data qubit is in the origin with total circle time of tau!
        # cOffset allows to move not perfectly above the data qubit
        r=circ_motion(t, args['cOpts'])
    return args['Hz'] + args['J']/norm(r)**3 * ( 1. -3. * r[2]**2/norm(r)**2 ) * args['Hd'] + args['J']/norm(r)**3 * ( 2. -3.*r[0]**2/norm(r)**2 -3.*r[1]**2/norm(r)**2 ) * (phase*args['H12'] + conjugate(phase)*args['H21'])



# Circular motion time dependence
# total time tau (for whole circle like in paper)
tau=3.3e-3 #e.g. 1ms

# GEOMETRY
# size of simulation square, D=ata qubit separation, d=sheet separation
D=400e-9
d=40e-9

# constants
ge=1.99875
gP=ge-2.5e-4
gBi=2.0003
ABi=1475e6
AP=118e6


muB=9.27e-24
Bfield=300e-3/cp.hbar   # division by hbar needed because not implemented in mesolve lindblad?
g1=gBi#10
g2=gP
J=cp.mu_0*ge**2 * muB**2/(4*pi)/cp.hbar
Delta=(g2-g1)*muB*Bfield

#
sigma2z = tensor(qeye(2), sigmaz())
Hz= Delta*sigma2z

# needed time independent diagonal and off diagonal parts. diagonal part is just sigmazz. off diagonal part is H_12 and H_21 part of the full matrix.
args={'Hz': Hz, 'Hd': tensor(sigmaz(), sigmaz()), 'H12': Qobj([[0,0,0,0],[0,0,0,0],[0,1,0,0],[0,0,0,0]]), 'H21': Qobj([[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]])}

# specify simulation arguments
args.update({'J': J, 'Delta': Delta, 'r': array([0,0,d]), 'circ': True, 'cOpts': {'pJit': False, 'rstd': 1.e-9, 'cOffset': zeros(3), 'tau': tau, 'd': d, 'D': D} })

# intial state
#               probe |+>               data
psi0 = tensor(qubit_state(pi/2.,0), qubit_state(0.,0,))  

# use time independent staying on top of each other 2*78e-6 does the pi/2 rotation we want!
#tlist=linspace(0, 2*78e-6, 100) #one simulation is only a quarter of the turn!
tlist=linspace(0, tau/4., 300)#40000) #one simulation is only a quarter of the turn!
result = mesolve(H_RWA, psi0, tlist, [], [], args)


#qsave(result.states, 'states')

# Plot Bloch
db=Bloch()

for t in range(0,len(result.states),10):
    db.add_states(result.states[t].ptrace(0), kind='point')
    db.add_states(result.states[t].ptrace(1), kind='point')
db.show()

# Test plot circ motion
"""
r=circ_motion(tlist, args['cOpts'])
plot(r[0], r[1])
show()
"""
