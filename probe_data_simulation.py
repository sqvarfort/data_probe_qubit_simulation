from H_RWA import H_RWA
from Plotter import Plotter

from DataQubitDisplacement import *
from Simulation import Simulation

from qutip import *


# here abrupt movement. tau is the time spend in interaction with each individual qubit
ham=H_RWA("ham_conf.yml")
tau=7.7155254093495203e-05

hamiltonian = ham.getHfunc()
initial_states = [(pi/2.0,0),(0,0),(0,0),(pi,0),(0,0)] # [(theta,phi),...] for all qubits
mesolve_args = ham.getArgs()

sim = Simulation(hamiltonian,initial_states,mesolve_args)

sim.run(4*tau, 4.*200)

result = sim.last_run_all

# Plot Bloch
db=Bloch()

for t in range(0,len(result),10):
    db.add_states(result[t].ptrace(0), kind='point')
    #db.add_states(result.states[t].ptrace(1), kind='point')
db.show()

