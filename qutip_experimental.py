from qutip import *
from pylab import *
from scipy import constants as cp
from numpy import random


def qubit_state(theta,phi):
    return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])

# define needed operators
sigma1z  = tensor(sigmaz(), qeye(2))
sigma2z = tensor(qeye(2), sigmaz())
sigmaxx = tensor(sigmax(), sigmax())
sigmayy = tensor(sigmay(), sigmay())
sigmazz = tensor(sigmaz(), sigmaz())
sigmaxy = tensor(sigmax(), sigmay())
sigmaxz = tensor(sigmax(), sigmaz())
sigmayx = tensor(sigmay(), sigmax())
sigmayz = tensor(sigmay(), sigmaz())
sigmazx = tensor(sigmaz(), sigmax())
sigmazy = tensor(sigmaz(), sigmay())


def circ_motion(t, args):
    r=args['cOffset'] + array( [ -args['D']/2. + args['D']/sqrt(2)*sin(2*pi/args['tau'] * t)  ,  -args['D']/2. + args['D']/sqrt(2)*cos(2*pi/args['tau'] * t)  ,  args['d'] ] )
    #allow random path jitter (simulates not perfect movement of mems stage) ~~1nm?
    if args['pJit']:
        # add gavins function here!
        r+=random.normal(0, args['rstd'], 3)
        # use below for test plot
        #r=args['cOffset'] + array( [ -args['D']/2. + random.normal(args['D']/sqrt(2)*sin(2*pi/args['tau'] * t), args['rstd'])  ,  -args['D']/2. + random.normal(args['D']/sqrt(2)*cos(2*pi/args['tau'] * t), args['rstd'])  ,  random.normal(args['d'], args['rstd']) ] )
    return r

# implements the whole time dependence
def H_full(t, args):
    r=args['r']
    if args['circ']:
        #implements circular motion from with radius D/sqrt2 where D is the separation between donors. data qubit is in the origin with total circle time of tau!
        # cOffset allows to move not perfectly above the data qubit
        r=circ_motion(t, args['cOpts'])
    return args['Hz'] + args['J']/norm(r)**3 * args['Hxyz'] -3*args['J']/norm(r)**3 * ( r[0]**2/norm(r)**2 * args['Hxx'] + r[1]**2/norm(r)**2 * args['Hyy'] + r[2]**2/norm(r)**2  * args['Hzz'] + r[0]*r[1]/norm(r)**2  * args['Hxy'] + r[0]*r[2]/norm(r)**2  * args['Hxz'] + r[1]*r[2]/norm(r)**2  * args['Hyz'] )



# Circular motion time dependence
# total time tau (for whole circle like in paper)
tau=1.2e-3 #e.g. 1ms

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

# go into probe Zeeman referemce frame, here proe is 1st qubit!
g=g1
g1=g1-g
g2=g2-g

# Time independent Hamilktonians. They need to be passed via args variable! becasue every qobj is converted to sparse matrix see mesolve man
Hz= muB * Bfield * ( g1*sigma1z + g2*sigma2z  )
Hxy=(sigmaxy + sigmayx)
Hxz=(sigmaxz + sigmazx)
Hyz=(sigmayz + sigmazy)
# add all needed time independent hamiltonians to args
args={'Hz': Hz, 'Hxyz': sigmaxx+sigmayy+sigmazz, 'Hxx': sigmaxx, 'Hyy': sigmayy, 'Hzz': sigmazz, 'Hxy': Hxy, 'Hxz': Hxz, 'Hyz': Hyz}

# specify simulation arguments
#args.update({'J': J, 'circ': False, 'pJit': False, 'cOffset': zeros(3), 'tau': tau, 'r': array([0,0,d]), 'd': d, 'D': D})
args.update({'J': J, 'Delta': (g2-g1)*muB*Bfield, 'r': array([0,0,d]), 'circ': False, 'cOpts': {'pJit': True, 'rstd': 1.e-9, 'cOffset': zeros(3), 'tau': tau, 'd': d, 'D': D} })

# intial state
#               probe |+>               data
psi0 = tensor(qubit_state(pi/2.,0), qubit_state(pi,0,))  

# use time independent staying on top of each other 2*78e-6 does the pi/2 rotation we want!
#tlist=linspace(0, 2*78e-6, 3000) #one simulation is only a quarter of the turn!
#result = mesolve(H_full(0, args), psi0, tlist, [], [])
#result = mesolve(H_full, psi0, tlist, [], [], args)


# use time dependent
tlist=linspace(0, tau/4., 100)#40000) #one simulation is only a quarter of the turn!
#result = mesolve(H, psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))

#qsave(result.states, 'states')

# Plot Bloch
"""
db=Bloch()

for t in range(0,len(result.states),100):
    db.add_states(result.states[t].ptrace(0), kind='point')
    db.add_states(result.states[t].ptrace(1), kind='point')
db.show()
"""

# Test plot circ motion
"""
r=circ_motion(tlist, args['cOpts'])
plot(r[0], r[1])
show()
"""
