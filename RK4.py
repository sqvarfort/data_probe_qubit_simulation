""" This is an RK4 solver that evolves a density matrix under some Hamiltonian, to be replaced by the full Lindblad equation.
Changes:
- Change density matrix to look at Bloch vector? Better intuition... But bad for future two-qubit states, can't have 4D Bloch vector...
- Implementing real Hamiltnonians: need to think of a way we can add error terms in a modular way. Suggestion: write long expression with coefficients that can be turned 'on' and 'off'. 
"""

import numpy as np
from numpy import exp


def RK4(Lindblad, rho, h, time): # RK4 solver.
    k1 = Lindblad(time, rho)
    k2 = Lindblad(time + h/2., rho + h/2.*k1)
    k3 = Lindblad(time + h/2., rho + h/2.*k2)
    k4 = Lindblad(time + h, rho + h*k3)
    return (k1 + 2.*k2 + 2.*k3 + k4)*h/6.



sigmax = np.matrix([[0,1], [1,0]])
sigmaz = np.matrix([[1,0],[0,-1]])
rho = np.matrix([[1,1],[1,1]])*1/2. # Start in |+> state


def commutator(A,B):
    return np.dot(A,B) - np.dot(B,A)

def Heisenberg(time, rho): # Very simple Hamiltonian
    return -1j*commutator(sigmaz,rho)




import matplotlib.pyplot as plt
plt.plot()


h = 0.1 # Time-step
time = 0.0 # Starting time

# Iterate over RK4-method
for t in range(0,100):
    dy = RK4(Heisenberg, rho, h, time) # Invke solver
    rho = rho + dy
    time = time + h
    print rho
    plt.scatter(time, rho[0,1]) # Plot off-diagonal matrix element
    plt.scatter(time, rho[1,0], color = 'red') # Plot the other off-diagonal element


plt.show()
