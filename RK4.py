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

#def class()


def RK4(Lindblad, rho, h, time, rData, rProbe): # RK4 solver.
    k1 = Lindblad(time, rho)
    k2 = Lindblad(time + h/2., rho + h/2.*k1)
    k3 = Lindblad(time + h/2., rho + h/2.*k2)
    k4 = Lindblad(time + h, rho + h*k3)
    return (k1 + 2.*k2 + 2.*k3 + k4)*h/6.



sigmax = np.matrix([[0,1], [1,0]])
sigmaz = np.matrix([[1,0],[0,-1]])
rho = np.matrix([[1,1],[1,1]])*1/2. # Start in |+> state

# Define the probe qubit as a function of distance r away from the data qubit. The angle theta sets the starting value

def initialise_data(x,y,z):
    return 0

def initialise_qubit(theta, phi, probe):
    return np.matrix([[cos[theta/2.], sin[theta]*exp[-1j*phi]],
                      [sin[theta]*exp[1j*phi]], sin[theta/2.]])



def commutator(A,B):
    return np.dot(A,B) - np.dot(B,A)

def Heisenberg(time, rho): # Very simple Hamiltonian
    return -1j*commutator(sigmaz,rho)

def Lindblad(time, rho, rProbe, rData):
    return 0



import matplotlib.pyplot as plt
plt.plot()


h = 0.01 # Time-step
time = 0.0 # Starting time
rProbe = 0
rData = 0

# Iterate over RK4-method
for t in range(0,100):
    dy = RK4(Heisenberg, rho, h, time, rProbe, rData) # Invke solver
    rho = rho + dy
    time = time + h
    print rho
    plt.scatter(time, rho[0,1]) # Plot off-diagonal matrix element
    plt.scatter(time, rho[1,0], color = 'red') # Plot the other off-diagonal element

# How would we move the probe about? One way is to use a 5-dimensional array for one single qubit. In the first two indices, we have the matrix. In the last three, we have the coordinates of its position. BUt how would we describe the entangled state? Maybe we could contain all the information about their closeness in the Hamiltonian? But how do we then use the solver?

# First, I must figure out where the time-dependence comes in here. I think that the h-parameter in the RK4-method seems to

# I'd like to associate r with the function rho, but I am not sure whether this is possible with the solver. Solution: we have to track the two positions separately, and then feed them as arguments into the Hamiltonian.

# We should probably calculate the error of the computation as well.


plt.show()
