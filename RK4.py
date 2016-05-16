""" This is an RK4 solver that evolves a density matrix under some Hamiltonian, to be replaced by the full Lindblad equation.
Changes:
- Change density matrix to look at Bloch vector? Better intuition... But bad for future two-qubit states, can't have 4D Bloch vector...
- Implementing real Hamiltnonians: need to think of a way we can add error terms in a modular way. Suggestion: write long expression with coefficients that can be turned 'on' and 'off'.

User:
- initial positions
- orbit type
- initial qubit
- time and time-steps
- turn Lindblad terms on and off
- control of d and D, d = distance between planes, D = lattice constant

Class:
-

"""

import numpy as np
from numpy import exp
from math import exp, sin, cos

rData = np.array([0,0,0])
rProbe = np.array([0,0,0])

sigmax = np.matrix([[0,1], [1,0]])
sigmaz = np.matrix([[1,0],[0,-1]])
rho = np.matrix([[1,1],[1,1]])*1/2. # Start in |+> state

# We don't need to use the lattice spacing. This is how closely the qubits are situated. The interesting thing will be the motion of the qubit.
# Must very carefully conisder the relationship between time and the time for a cycle
# The cycle time is the maximum time that we require


# Lattice spaceing of data qubits in mu m.

# cycTime / h = number of iterations

class data_probe_sim(cycTime, h, distance, separation, size, Bfield, g1, g2):
    self.rProbe = rProbe # 3-dim position vector initialised by user
    self.cycleTime = cycleTime
    self.h = h # time-step
    self.separation = separation # Separation between the overhead and the body, initial value
    self.distance = distance # Distance between probe and data qubit
    self.size = 1000 # This decides the mesh size. This would be one point per nm.
    it_no = cycTime/(h*4) # Number of iterations per cycle. We only execute 1/4 cycles.
    self.g1 = g1
    self.g2 = g2
    #self.velocity = it_no / new distance point
    # Distance will always be the same
    #Geometry
    rData = np.array([0, 0, 0]) # Position data qubit in the middle of the sample, which is also the origin. Need to keep this to model displacement of the data qubit from its intended position.
    rProbe = np.array([0, size/2., separation])
    distance = np.linalg.norm(rProbe - rData)
    rtheta = np.arccos(rProbe[2]/distance)
    rphi = np.arctan(rProbe[1]/rProbe[0])

    self.Bfield = Bfield # Including the Bohr magneton
    muB = 1 # Set to unity, Only relationship that matters is relationship between J and B
    J = 1

    sigmax = np.matrix([[0,1], [1,0]])
    sigmay = np.matrix([[0,-1j], [1j, 0]])
    sigmaz = np.matrix([[1,0],[0,-1]])
    identity = np.matrix([[1,0],[0,1]])
    sigmavec = np.array([sigmax, sigmay, sigmaz])

    # Define the Hamiltonian


    def update_circ_orbit(): # This takes the current position rProbe and updates the variable
        rProbe = rProbe + 0 # Have fun....
        distance = np.linalg.norm(rProbe - rData)
        rtheta = np.arccos(rProbe[2]/distance)
        rphi = np.arctan(rProbe[1]/rProbe[0])

    def circ_obrit(rProbe): # Increment displacement
        return 0

    def update_geometry():
        distance = np.linalg.norm(rProbe - rData)
        rtheta = np.arccos(rProbe[2]/distance)
        rphi = np.arctan(rProbe[1]/rProbe[0])


    def lin_orbit(): # jumpy orbit displacement
        return 0

    def Hamiltonian(): # Make sure to update the distance before H is calculated at each point
        return muB * Bfield (g1 * tensordot[sigmaz, identity] + g2*tensordot[identity, sigmaz]) + J/(power[distance, 3])*(tensordot[sigmax, sigmax] + tensordot[sigmay, sigmay] + tensordot[sigmaz, sigmaz])
        - 3*(sin[rtheta]*cos[rphi]*tensordot[sigmax, identity] + sin[rtheta]*sin[rphi]*tensordot[sigmay, identity] + cos[rtheta]*tensordot[sigmaz, identity])*(sin[rtheta]*cos[rphi]*tensordot[identity, sigmax] + sin[rtheta]*sin[rphi]*tensordot[identity, sigmay] + cos[rtheta]*tensordot[identity, sigmaz])

# Removed time-dependence, could be that we need to put it back in
    def Lindblad(system):
        return  - 1j*commutator[Hamiltonian(), system]


    def RK4(Lindblad): # RK4 solver. Get global variables from class
        k1 = Lindblad(system)
        k2 = Lindblad(rho + h/2.*k1)
        k3 = Lindblad(rho + h/2.*k2)
        k4 = Lindblad(rho + h*k3)
        return (k1 + 2.*k2 + 2.*k3 + k4)*h/6.




# Define the probe qubit as a function of distance r away from the data qubit. The angle theta sets the starting value

def init_system():
        return tensordot(pRho, dRho)

def init_qubit(theta, phi):
    self.pRho = np.matrix([[cos[theta/2.], sin[theta]*exp[-1j*phi]],
                      [sin[theta]*exp[1j*phi]], sin[theta/2.]])



def commutator(A,B):
    return np.dot(A,B) - np.dot(B,A)

def Heisenberg(time, rho): # Very simple Hamiltonian
    return -1j*commutator(sigmaz,rho)


def Lindblad():

    # First implement the simple Hamiltonian

    return 0



import matplotlib.pyplot as plt
plt.plot()


h = 0.01 # Time-step
rProbe = 0
rData = 0

# Iterate over RK4-method
    def iterate():
        for t in range(0,it_no):
            dy = RK4(Heisenberg) # Invke solver
            system = system + dy # This should also be a global variables
            time = time + h #Only useful for the plotting
            rProbe = circ_orbi(rProbe)
            print system
            plt.scatter(time, system[0,1]) # Plot off-diagonal matrix element
            plt.scatter(time, system[1,0], color = 'red') # Plot the other off-diagonal element
            rProbe = update_circ_orbit(rProbe)
            update_geometry()


# How would we move the probe about? One way is to use a 5-dimensional array for one single qubit. In the first two indices, we have the matrix. In the last three, we have the coordinates of its position. BUt how would we describe the entangled state? Maybe we could contain all the information about their closeness in the Hamiltonian? But how do we then use the solver?

# First, I must figure out where the time-dependence comes in here. I think that the h-parameter in the RK4-method seems to

# I'd like to associate r with the function rho, but I am not sure whether this is possible with the solver. Solution: we have to track the two positions separately, and then feed them as arguments into the Hamiltonian.

# We should probably calculate the error of the computation as well.


plt.show()
