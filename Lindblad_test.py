import numpy as np
from pylab import *
from scipy import constants as cp
from qutip import *

from Lindblad import Lindblad

def qubit_state(theta,phi):
    return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])

# simple 1-qubit system
#rho0 = qubit_state(pi/2,pi) # theta,phi

# simple 2-qubit system
rho0 = tensor(qubit_state(pi/2,0), qubit_state(pi/2,0,)) 

# Hamiltonian
H = 1*qeye(4)

samples = np.linspace(0.0,100,50)

lind = Lindblad(2) # Collapse operators over (x) qubits
lind.dephasing(0.05,0) # Dephasing on first qubit (rate, qubit)
linds = lind.lindblads


#result = mesolve(H, rho0, samples, [], [])
result = mesolve(H, rho0, samples, linds, [])

plt.plot(samples, 0.5-0.5*(expect(tensor(sigmax(),qeye(2)),result.states)))
plt.show()


db=Bloch()

for t in range(len(result.states)):
    db.add_states(result.states[t].ptrace(0), kind='point')
    #db.add_states(result.states[t].ptrace(1), kind='point') # Current simulation only one qubit
db.show()