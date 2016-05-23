from qutip import *
from pylab import *
from scipy import constants as cp
from numpy import random


def qubit_state(theta,phi):
    return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])

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
print J


# GEOMETRY
# size of simulation square
size=400e-9
separation=40e-9
# vector between the spins   NEEDS TO BE NORMALISED i nHamiltonian!
rOffset=array([-size/2.,-size/2, separation])
r=array([0,0,separation])

# intial state
#               probe |+>               data
psi0 = tensor(qubit_state(pi/2.,0), qubit_state(pi,0,))
# operators
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


# Circular motion time dependence
# total time tau (for whole circle like in paper)
tau=1.2e-3 #e.g. 1ms

# Standard deviation for path fluctuations
rstd = 1e-9

def calc_r(t):
    r=rOffset+array([size/sqrt(2)*sin(2*pi/tau * t), size/sqrt(2)*cos(2*pi/tau * t), 0.])
    r_rand = rOffset + array([random.normal(size/sqrt(2)*sin(2*pi/tau * t), rstd), random.normal(size/sqrt(2)*cos(2*pi/tau * t), rstd), random.normal(0,rstd)])
    return r_rand

test_time = 0

for i in range(0,100):
    calc_r(test_time)


def Hxyz_coeff(t, args):
    r = calc_r(t)
    return J/norm(r)**3

Hxyz=sigmaxx + sigmayy + sigmazz

def Hxx_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[0]**2

Hxx=sigmaxx

def Hyy_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[1]**2

Hyy=sigmayy

def Hzz_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[2]**2

Hzz=sigmazz

def Hxy_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[0]*r[1]

Hxy=sigmaxy

def Hxz_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[0]*r[2]

Hxz=sigmaxz

def Hyx_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[0]*r[1]

Hyx=sigmayx

def Hyz_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[1]*r[2]

Hyz=sigmayz

def Hzx_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[0]*r[2]

Hzx=sigmazx

def Hzy_coeff(t, args):
    r=calc_r(t)
    return -3*J/norm(r)**5 * r[2]*r[1]

Hzy=sigmazy

# go into probe Zeeman referemce frame, here proe is 1st qubit!
g=g1
g1=g1-g
g2=g2-g

# Zeeman Hamitonian
Hz= muB * Bfield * ( g1*sigma1z + g2*sigma2z  )

Hi= J/norm(r)**3 * ( sigmaxx + sigmayy + sigmazz -3/norm(r)**2 * ( r[0]**2 * sigmaxx + r[1]**2 * sigmayy + r[2]**2 * sigmazz + r[0]*r[1] * sigmaxy + r[0]*r[2] * sigmaxz + r[1]*r[0]* sigmayx + r[1]*r[2] * sigmayz + r[2]*r[0] * sigmazx + r[2]*r[1]* sigmazy  ) )


# Full Hamiltonian time dependent
H=[Hz,[Hxyz, Hxyz_coeff],[Hxx, Hxx_coeff],[Hyy, Hyy_coeff],[Hzz, Hzz_coeff],[Hxy, Hxy_coeff],[Hxz, Hxz_coeff],[Hyx, Hyx_coeff],[Hyz, Hyz_coeff],[Hzx, Hzx_coeff],[Hzy, Hzy_coeff]]


#Htest=tensor(qeye(2), sigmaz())


# use time independent staying on top of each other 2*78e-6 does the pi/2 rotation we want!
tlist=linspace(0, 2*78e-6, 3000) #one simulation is only a quarter of the turn!
result = mesolve(Hz+Hi, psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))

# use time dependent
#tlist=linspace(0, tau/4., 40000) #one simulation is only a quarter of the turn!
#result = mesolve(H, psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))

qsave(result.states, 'states')

db=Bloch()

for t in range(0,len(result.states),100):
    db.add_states(result.states[t].ptrace(0), kind='point')
    db.add_states(result.states[t].ptrace(1), kind='point')
db.show()
