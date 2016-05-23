from qutip import *
import numpy as np
from pylab import *
from scipy import constants as cp


def qubit_state(theta,phi):
    return Qobj([[np.cos(theta/2.)], [np.exp(1j*phi) * np.sin(theta/2.)]])


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

def H_circ(t, args):
    #r=rOffset+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    #r=np.array([0,0,40e-9])
    #r=args['r']
    #Hxyz=(args['xx'] + args['yy'] + args['zz'])
    #Hxx=r[0]**2 * args['xx']
    #Hyy=r[1]**2 * args['yy']
    #Hzz=r[2]**2 * args['zz']
    #Hxy=r[0]*r[1] * (args['xy'] + args['yx'])
    #Hxz=r[0]*r[2] * (args['xz'] + args['zx'])
    #Hyz=r[1]*r[2] * (args['yz'] + args['zy'])
    Hxyz=(args['xx'] + args['yy'] + args['zz'])
    #Hxx=r[0]**2 * args['xx']
    #Hyy=r[1]**2 * args['yy']
    Hzz=args['sep']**2 * args['zz']
    #Hxy=r[0]*r[1] * (args['xy'] + args['yx'])
    #Hxz=r[0]*r[2] * (args['xz'] + args['zx'])
    #Hyz=r[1]*r[2] * (args['yz'] + args['zy'])
    return args['Hz'] + args['J']/args['sep']**3 * Hxyz -3*args['J']/args['sep']**5 * ( Hzz )
    #return args['Hz'] + args['J']/np.linalg.norm(r)**3 * Hxyz -3*args['J']/np.linalg.norm(r)**5 * ( Hxx + Hyy + Hzz + Hxy + Hxz + Hyz )

def Hxyz_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    return args['J']/norm(r)**3 

Hxyz=sigmaxx + sigmayy + sigmazz 

def Hxx_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    return -3*J/norm(r)**5 * r[0]**2 + Hxyz_coeff(t,args)

Hxx=sigmaxx

def Hyy_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    return -3*J/norm(r)**5 * r[1]**2 + Hxyz_coeff(t,args)

Hyy=sigmayy

def Hzz_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    return -3*J/norm(r)**5 * r[2]**2 + Hxyz_coeff(t,args)

Hzz=sigmazz

def Hxy_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    r=rOffset+array([size/sqrt(2)*sin(2*pi/tau * t), size/sqrt(2)*cos(2*pi/tau * t), 0.])
    return -3*J/norm(r)**5 * r[0]*r[1]

Hxy=sigmaxy
Hyx=sigmayx

def Hxz_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    r=rOffset+array([size/sqrt(2)*sin(2*pi/tau * t), size/sqrt(2)*cos(2*pi/tau * t), 0.])
    return -3*J/norm(r)**5 * r[0]*r[2]

Hxz=sigmaxz
Hzx=sigmazx

def Hyz_coeff(t, args):
    r=args['r']
    if args['circ']:
        r=args['rStart']+array([size/sqrt(2)*sin(2*pi/args['tau'] * t), size/sqrt(2)*cos(2*pi/args['tau'] * t), 0.])
    return -3*J/norm(r)**5 * r[1]*r[2]

Hyz=sigmayz
Hzy=sigmazy

def ttest(t,args):
    return 1.
def H_test(t, args):
    #return 1.*args['Hz']+1.*args['Hi']
    return args['Hz'] + args['A'] * (args['xyz'] ) + args['B'] * ( args['zz'] )

#-----------------------------------------

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
J=cp.mu_0*ge**2 * muB**2/(4*np.pi)/cp.hbar

# GEOMETRY
# size of simulation square
size=400e-9
separation=40e-9
# vector between the spins   NEEDS TO BE NORMALISED i nHamiltonian!
rOffset=np.array([-size/2.,-size/2, separation])
r=np.array([0,0,separation])

# Circular motion time dependence
# total time tau (for whole circle like in paper)
tau=1.2e-3 #e.g. 1ms


# go into probe Zeeman referemce frame, here proe is 1st qubit!
g=g1
g1=g1-g
g2=g2-g

# Zeeman Hamitonian
Hz= muB * Bfield * ( g1*sigma1z + g2*sigma2z  )

Hi= J/norm(r)**3 * ( sigmaxx + sigmayy + sigmazz -3/norm(r)**2 * ( r[0]**2 * sigmaxx + r[1]**2 * sigmayy + r[2]**2 * sigmazz + r[0]*r[1] * sigmaxy + r[0]*r[2] * sigmaxz + r[1]*r[0]* sigmayx + r[1]*r[2] * sigmayz + r[2]*r[0] * sigmazx + r[2]*r[1]* sigmazy  ) )


# Full Hamiltonian time dependent
#H=[Hz,[Hxx, Hxx_coeff],[Hyy, Hyy_coeff],[Hzz, Hzz_coeff],[Hxy, Hxy_coeff],[Hyx, Hxy_coeff],[Hxz, Hxz_coeff],[Hzx, Hxz_coeff],[Hyz, Hyz_coeff], [Hzy, Hyz_coeff]]

H=[Hz, [Hxyz, 'A'], [Hzz, 'B']]

# intial state
#               probe |+>               data
psi0 = tensor(qubit_state(np.pi/2.,0), qubit_state(0,0,))    

print {'Hz': Hz, 'Hi': Hi, 'xyz': Hxyz, 'zz': Hzz, 'A': J/separation**3, 'B': -3.*J/separation**5}
# use time independent staying on top of each other 2*78e-6 does the pi/2 rotation we want!
tlist=np.linspace(0, 2*78e-6, 3000) #one simulation is only a quarter of the turn!
#result = mesolve(Hz+Hi, psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))
#result = mesolve(H_test, psi0, tlist, [], [], {'Hz': Hz, 'Hi': Hi, 'xyz': Hxyz, 'zz': Hzz, 'A': J/separation**3, 'B': -3.*J/separation**5})#, options=Odeoptions(nsteps=100000))
#result = mesolve([Hz,[Hi,ttest]], psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))
args={'xx': sigmaxx, 'yy': sigmayy, 'zz': sigmazz, 'xy': sigmaxy, 'yx': sigmayx, 'zx': sigmazx, 'yz': sigmayz, 'xz': sigmaxz, 'zy': sigmazy, 'Hz': Hz, 'J': J, 'circ': False, 'tau': tau, 'rStart': rOffset, 'r': np.array([0,0,separation]), 'sep': 40e-9}
#Hc=H_circ(0,args)
#result = mesolve(Hc, psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))
#result = mesolve(H_circ, psi0, tlist, [], [],args)#, options=Odeoptions(nsteps=100000))
#result = mesolve(H, psi0, tlist, [], [], {'A': J/separation**3, 'B': -3.*J/separation**5})#, options=Odeoptions(nsteps=100000))

# use time dependent
#tlist=linspace(0, tau/4., 40000) #one simulation is only a quarter of the turn!
#result = mesolve(H, psi0, tlist, [], [])#, options=Odeoptions(nsteps=100000))

#qsave(result.states, 'abrupt-0')
"""
db=Bloch()

for t in range(0,len(result.states),100):
    db.add_states(result.states[t].ptrace(0), kind='point')
    db.add_states(result.states[t].ptrace(1), kind='point')
db.show()
""" 
