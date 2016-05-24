from H_RWA import H_RWA
from Plotter import Plotter

from DataQubitDisplacement import *
from Simulation import Simulation

from qutip import *


# here abrupt movement. tau is the time spend in interaction with each individual qubit
ham=H_RWA("ham_conf.yml")
tau=7.7155254093495203e-05
info = 'something'

hamiltonian = ham.getHfunc()
initial_states = [(pi/2.0,0),(0,0),(0,0),(pi,0),(0,0)] # [(theta,phi),...] for all qubits
mesolve_args = ham.getArgs()
iterations = 200

sim = Simulation(hamiltonian,initial_states,mesolve_args)

sim.run(tau, iterations, 1)

result = sim.last_run_all

# Plot Bloch
db=Bloch()
colors = ["g"]#,"r","g","#CC6600"]#,"r","g","#CC6600"]
db.point_color = colors
db.point_marker = ['o']

for t in range(0,len(result),10):
    db.add_states(result[t].ptrace(0), kind='point')
    #db.add_states(result.states[t].ptrace(1), kind='point')
db.show()

Plotter(result, info)


"""

for i in range(0,2):
    sim.run(tau,iterations)
    result_states = sim.last_run_all
    step_data = sim.last_run_quarter_cycle
    final_states.append(sim.last_run_all[-1])
    sim.reset_system_state()

qsave(result_states, 'states')
"""
