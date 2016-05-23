from qutip import *
from pylab import *
from scipy import constants as cp

from Lindblad import Lindblad

class Simulation():
    def __init__(self,hamiltonian,states,mesolve_args):
        self.hamiltonian = hamiltonian
        self.args = mesolve_args

        self.lind = Lindblad(2) # Add Lindblads with eg Simulation.lind.dephasing(1.0e3,0)
        
        self.initial_states = states
        self.full_state = self.set_system_state(states)
        
        #self.start_states = [] # Stores the initial 5-qubit Qobjs before every run
        self.final_states = [] # Stores the final 5-qubit Qobjs after every run
        self.last_run_all = []
        self.last_run_quarter_cycle = []
        self.last_run_metadata = self._store_metadata(None,None,None)  
        
    def qubit_state(self,theta,phi):
        return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])
    
    def set_system_state(self, states):
        # states should be of form [(theta,phi),...] for all qubits
        q_probe = self.qubit_state(states[0][0],states[0][1])
        q_data1 = self.qubit_state(states[1][0],states[1][1])
        q_data2 = self.qubit_state(states[2][0],states[2][1])
        q_data3 = self.qubit_state(states[3][0],states[3][1])
        q_data4 = self.qubit_state(states[4][0],states[4][1])
        return tensor(q_probe,q_data1,q_data2,q_data3,q_data4)
    
    def _run_quarter_cycle(self,state,time,steps):
        tlist=linspace(0, time, steps) #one simulation is only a quarter of the turn!
        it_res = mesolve(self.hamiltonian, state, tlist, self.lind.lindblads, self.args)
        for state in it_res.states: self.last_run_all.append(state)
        return it_res.states[-1] # return last Qobj
        
    def _store_metadata(self,time,steps,cycles):
        self.last_run_metadata = {'cycles':cycles, 'time':time, 'steps':steps, 'lindblad':self.lind.lindblads, 'initial_states':self.initial_states, 'args':self.args}
        return self.last_run_metadata
        
    def run(self,time,steps,cycles=4):
        self.last_run_quarter_cycle = []
        self.last_run_all = []
        cycle_time = time/cycles
        cycle_steps = steps/cycles
        
        #self.initial_states.append(self.full_state) # Optionally store initial states
        
        for cycle in range(cycles):
            probe_qubit = self.full_state.ptrace(0)
            data_qubit = self.full_state.ptrace(cycle+1)
            
            state = tensor(probe_qubit,data_qubit)
            state = self._run_quarter_cycle(state,cycle_time,cycle_steps)
            self.last_run_quarter_cycle.append(state)
            
            all_data_qubits = [self.full_state.ptrace(i+1) for i in range(4)]
            probe_qubit = state.ptrace(0)
            all_data_qubits[cycle] = state.ptrace(1)
            
            self.full_state = tensor(probe_qubit,all_data_qubits[0],all_data_qubits[1],all_data_qubits[2],all_data_qubits[3])
        
        self.final_states.append(self.full_state)
        self._store_metadata(time,steps,cycles)
        
        return self.full_state