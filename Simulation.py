from qutip import *
from pylab import *
from scipy import constants as cp

from Lindblad import Lindblad

class Simulation():
    def __init__(self,hamiltonian):
        self.hamiltonian = hamiltonian
        self.final_states = []
        self.last_run = []
        self.lind = Lindblad(2) # Add Lindblads with eg Simulation.lind.dephasing(1.0e3,0)
        
        initial_states = [(pi/2.0,0),(0,0),(0,0),(0,0),(0,0)]
        self.full_state = self.full_system_state(initial_states)
    
    def qubit_state(self,theta,phi):
        return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])
    
    def full_system_state(self, states):
        # states should be of form [(theta,phi),...] for all qubits
        q_probe = self.qubit_state(states[0][0],states[0][1])
        q_data1 = self.qubit_state(states[1][0],states[1][1])
        q_data2 = self.qubit_state(states[2][0],states[2][1])
        q_data3 = self.qubit_state(states[3][0],states[3][1])
        q_data4 = self.qubit_state(states[4][0],states[4][1])
        return tensor(q_probe,q_data1,q_data2,q_data3,q_data4)
    
    def _run_quarter_cycle(self,state,time=2*78e-6,steps=3000):
        tlist=linspace(0, time, steps) #one simulation is only a quarter of the turn!
        it_res = mesolve(self.hamiltonian, state, tlist, self.lind.lindblads, [])
        for state in it_res.states: self.last_run.append(state)
        return it_res.states[-1] # return last Qobj
        
    def run_cycles(self,cycles,state,time=8*78e-6,steps=12000):
        self.last_run = []
        cycle_time = time/cycles
        cycle_steps = steps/cycles
        
        
        
        for cycle in range(cycles):
            
            state = self._run_quarter_cycle(state,cycle_time,cycle_steps)
        self.final_states.append(state)
        return state