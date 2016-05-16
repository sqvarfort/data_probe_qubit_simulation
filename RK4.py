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
from math import sin, cos, sqrt, pi
from cmath import exp
import matplotlib.pyplot as plt
import qutip

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

class simulation(object):
    def __init__(self, cycleTime, iterations, separation, size):
        self.cycleTime = cycleTime
        self.iterations = iterations
        self.h = self.cycleTime / self.iterations
        self.separation = separation # Separation between the overhead and the body, initial value
        self.size = size # This decides the mesh size. This would be one point per nm.
        self.orbit_radius = self.size/sqrt(2.) # Calculate the orbit radius
        it_no = self.cycleTime/(self.h*4.) # Number of iterations per cycle. We only execute 1/4 cycles.

        self.g1,self.g2 = 1.,1.
        #self.velocity = it_no / new distance point
        # Distance will always be the same
        #Geometry
        self.rData = np.array([0, 0, 0]) # Position data qubit in the middle of the sample, which is also the origin. Need to keep this to model displacement of the data qubit from its intended position.
        self.rProbe = np.array([- self.size/2., self.orbit_radius - self.size/2., self.separation]) # PLaces
        self.distance = np.linalg.norm(self.rProbe - self.rData)
        self.rtheta = np.arccos(self.rProbe[2]/self.distance) # Work out angles between data-qubit and probe-qubit
        self.rphi = np.arctan(self.rProbe[1]/self.rProbe[0])

        self.Bfield =  1.0
        self.muB = 1. # Set to unity, Only relationship that matters is relationship between J and B
        self.J = 1.

        self.sigmax = np.matrix([[0,1], [1,0]])
        self.sigmay = np.matrix([[0,-1j], [1j, 0]])
        self.sigmaz = np.matrix([[1,0],[0,-1]])
        self.identity = np.matrix([[1,0],[0,1]])


        # Initialise qubits, all in the |0> state
        self.ptheta = pi/2.
        self.pphi = 0.
        self.dtheta = 0.0
        self.dphi = 0.0



    def update_circ_orbit(self, time): # Based on the time in the iteration
        #Radius will be at (-size/2, -size/2)
        self.rProbe = np.array([-self.size/2. + self.orbit_radius*sin(time*pi/2), -self.size/2. + orbit_radius*cos(time*pi/2), self.separation])

        #Something like this. This is not correct, need the correct starting point, but will follow a similar quarter-circle parametrisation.


    def update_geometry(self, rProbe, rData):
        self.distance = np.linalg.norm(rProbe - rData)
        self.rtheta = np.arccos(rProbe[2]/self.distance)
        self.rphi = np.arctan(rProbe[1]/rProbe[0])


    def Hamiltonian(self, muB, Bfield, g1, g2, J): # Make sure to update the distance before H is calculated at each point
        return muB * Bfield *(g1 * np.kron(self.sigmaz, self.identity) + g2*np.kron(self.identity, self.sigmaz))  + J/(self.distance**3) *        (np.kron(self.sigmax, self.sigmax) + np.kron(self.sigmay, self.sigmay) + np.kron(self.sigmaz, self.sigmaz) - 3*(sin(self.rtheta)*cos(self.rphi)*np.kron(self.sigmax, self.identity) + sin(self.rtheta)*sin(self.rphi)*np.kron(self.sigmay, self.identity) + cos(self.rtheta) * np.kron(self.sigmaz, self.identity))          *     (sin(self.rtheta)*cos(self.rphi) * np.kron(self.identity, self.sigmax) + sin(self.rtheta)*sin(self.rphi)*np.kron(self.identity, self.sigmay) + cos(self.rtheta)*np.kron(self.identity, self.sigmaz)))



    def RK4(self, Lindblad, system, h): # RK4 solver. Get global variables from class
        k1 = Lindblad(system)
        k2 = Lindblad(system + h/2.*k1)
        k3 = Lindblad(system + h/2.*k2)
        k4 = Lindblad(system + h*k3)
        return (k1 + 2.*k2 + 2.*k3 + k4)*h/6.




    # Define the probe qubit as a function of distance r away from the data qubit. The angle theta sets the starting value

    def initialise_qubit(self, theta, phi):
        return np.matrix([[cos(theta/2.)**2, (1/2.)*sin(theta)*exp(-1.0j*phi)],[(1/2.)*sin(theta)*exp(1.0j*phi), sin(theta/2.)**2]])

    def commutator(self, A,B):
        return np.dot(A,B) - np.dot(B,A)

    def Heisenberg(self, system): # Very simple Hamiltonian
        return -1j*self.commutator(np.kron(self.sigmaz, self.sigmaz),system)

    def Lindblad(self, system):
        return -1j*self.commutator(self.Hamiltonian(self.muB, self.Bfield, self.g1, self.g2, self.J), system)

    # First implement the simple Hamiltonian




    # Iterate over RK4-method
    def evolve_system(self):
        time = 0.0
        ProbeQubit = self.initialise_qubit(self.ptheta, self.pphi) # Initialise qubits
        print ProbeQubit
        DataQubit = self.initialise_qubit(self.dtheta, self.dphi)
        system = np.kron(ProbeQubit, DataQubit) #Starting state - separable
        plt.plot()
        # Start iteration
        for i in range(0,int(self.iterations)):
            dy = self.RK4(self.Lindblad, system, self.h) # Invke solver
            system = system + dy # This should also be a global variables
            time = time + self.h #Only useful for the plotting
            plt.scatter(time, system[0,2], color = 'blue') # Plot the other off-diagonal element
            plt.scatter(time, system[2,0], color = 'red')
            #self.update_circ_orbit(time) # Calculate new point in orbit
            #self.update_geometry(rProbe, rData) # update the other geomtry parameters
            # Ideally, we only want the Hamiltonian and other Lindblad operators to be calculated once.
        # Trace out system and store it
        plt.show()

    def plotting(self):
        return 0

    def partial_trace(self):
        return 0


    # How would we move the probe about? One way is to use a 5-dimensional array for one single qubit. In the first two indices, we have the matrix. In the last three, we have the coordinates of its position. BUt how would we describe the entangled state? Maybe we could contain all the information about their closeness in the Hamiltonian? But how do we then use the solver?

    # First, I must figure out where the time-dependence comes in here. I think that the h-parameter in the RK4-method seems to

    # I'd like to associate r with the function rho, but I am not sure whether this is possible with the solver. Solution: we have to track the two positions separately, and then feed them as arguments into the Hamiltonian.

    # We should probably calculate the error of the computation as well.



    def run_full_cycle(): #Runs four iterations.
        for i in range(0,3):
            system = tensordot(ProbeQubit, DataQubit)
            evolve_system()
            # Include plotting stuff
            # Trace out probe qubit
            DataQubit = initialise_qubit(dtheta, dphi)


# Maybe it would make more sense to turn this into two classes - one class that creates and handles the system, and one that performs operations on it.

# try to enter units in nm
class_object = simulation(5.0, 100.0, 40, 10000)

class_object.evolve_system()
