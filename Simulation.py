import numpy as np
from qutip import *
from pylab import *
from scipy import constants as cp


from Lindblad import Lindblad
from DataQubitDisplacement import *

class Simulation():
    '''
    Prepare and run simulations on 5-qubit system using given hamiltonian
    
    Run simulation with Simulation.run(time,steps)
    
    This class stores some data after each run, some of which is
    persistent over several run()s, some of which gets
    overwritten after each run().
    
    Stored data:
    - final_states (each run() simulation appends its final 5-qubit state here) [PERSISTENT]
    - last_run_all (stores each 2-qubit state after each step of the simulaiton) [OVERWRITTEN]
    - last_run_quarter_cycle (stores the 2-qubit state after simulation of each quandarnt) [OVERWRITTEN]
    - last_run_metadata (stores some metadata from the last run) [OVERWRITTEN]
    '''
    def __init__(self,hamiltonian,states,mesolve_args):
        '''
        Prepare a simulation with given hamiltonian, mesolve arguments, and initial states
        '''
        self.hamiltonian = hamiltonian
        self.args = mesolve_args

        self.lind = Lindblad(2) # Add Lindblads with eg Simulation.lind.dephasing(1.0e3,0)
        
        self.initial_states = states
        self.set_system_state(states)
        
        self.progress_bar = True
        
        self.qubit_offsets_bool = False
        self.qubit_offsets_list = [[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0]]
        self.qubit_offsets_settings = {'radius':0.0, 'half_height':0.0, 'type':'uniform', 'convert':False}
        self.add_qubit_offsets = True
        
        twirl1 = tensor(qeye(2),qeye(2),qeye(2),qeye(2),qeye(2)) # IIII
        twirl2 = tensor(qeye(2),qeye(2),qeye(2),sigmax(),sigmax()) # IIXX
        twirl3 = tensor(qeye(2),qeye(2),sigmax(),qeye(2),sigmax()) # IXIX
        twirl4 = tensor(qeye(2),qeye(2),sigmax(),sigmax(),qeye(2)) # IXXI
        self.twirls = [twirl1,twirl2,twirl3,twirl4]
        self.twirls_str = ['IIII','IIXX','IXIX','IXXI']
        self.twirl = None
        
        #self.start_states = [] # Stores the initial 5-qubit Qobjs before every run
        self.final_states = [] # Stores the final 5-qubit Qobjs after every run
        self.last_run_all = []
        self.last_run_quarter_cycle = []
        self.last_run_metadata = self._store_metadata(None,None,None) 

        
        
    def qubit_state(self,theta,phi):
        return Qobj([[cos(theta/2.)], [exp(1j*phi) * sin(theta/2.)]])
    
    def set_system_state(self, states):
        '''
        Set 5-qubit system state using given rotations from |0>
        
        Probe qubit is qubit 0
        '''
        # states should be of form [(theta,phi),...] for all qubits
        q_probe = self.qubit_state(states[0][0],states[0][1])
        q_data1 = self.qubit_state(states[1][0],states[1][1])
        q_data2 = self.qubit_state(states[2][0],states[2][1])
        q_data3 = self.qubit_state(states[3][0],states[3][1])
        q_data4 = self.qubit_state(states[4][0],states[4][1])
        self.full_state = ket2dm( tensor(q_probe,q_data1,q_data2,q_data3,q_data4) )
        return self.full_state
        
    def reset_system_state(self):
        '''
        Resets the system back to its initial state
        '''
        self.set_system_state(self.initial_states)
        
    
    def _run_quarter_cycle(self,state,time,steps):
        tlist=linspace(0, time, steps) #one simulation is only a quarter of the turn!
        it_res = mesolve(self.hamiltonian, state, tlist, self.lind.lindblads, e_ops=[], args=self.args, progress_bar=self.progress_bar)
        for state in it_res.states: self.last_run_all.append(state)
        return it_res.states[-1] # return last Qobj
        
    def _store_metadata(self,time,steps,cycles):
        self.last_run_metadata = {'cycles':cycles, 'time':time, 'steps':steps, 'lindblad':self.lind.lindblads, 'initial_states':self.initial_states, 'args':self.args, 'qubit_offsets_list':self.qubit_offsets_list, 'twirl':self.twirl}
        return self.last_run_metadata
        
    def set_data_qubit_offsets(self,offset_list):
        '''
        Set a given list of 4 [x,y,z] qubit offsets for use in the simulation
        
        Clear offsets by running this function with the arguments None or False
        '''
        if offset_list: 
            self.qubit_offsets_list = offset_list
            self.qubit_offsets_bool = True
        else: # Allows passing None or False to turn off offsets
            self.qubit_offsets_bool = False
            self.qubit_offsets_list = [[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0]]
            
    def generate_data_qubit_offsets(self,radius,half_height,type='uniform',convert=False):
        '''
        Generate a list of 4 [x,y,z] qubit offsets representing
        random data qubit displacement
        
        Stores the settings used to generate offsets for use in 
        '''
        offset_list = []
        if type == 'uniform':
            for i in range(4):
                offset_list.append( random_uniform_cylinder(radius,half_height) )
        elif type == 'gaussian': # if using gaussian must specify standard deviations
            for i in range(4):   # alternatively convert uniform to sd with convert=True
                offset_list.append( random_gaussian_cylinder(radius,half_height,convert) )
        
        self.qubit_offsets_settings = {'radius':radius, 'half_height':half_height, 'type':type, 'convert':convert}
        
        self.qubit_offsets_list = offset_list
        self.qubit_offsets_bool = True
        
        return offset_list
        
    def regenerate_data_qubit_offsets(self):
        self.generate_data_qubit_offsets(self.qubit_offsets_settings['radius'],
                                         self.qubit_offsets_settings['half_height'],
                                         self.qubit_offsets_settings['type'],
                                         self.qubit_offsets_settings['convert'])
    
    def choose_twirl(self,i='default'):
        if not i: self.twirl = None # if call choose_twirl with False or None, apply no twirl
        elif type(i)==int: self.twirl = i # address twirls with 1,2,3,4; address array with 0,1,2,3
        elif i=='default': self.twirl = np.random.randint(4)+1 # No i specified, generate
        elif type(i)==bool and i: self.twirl = np.random.randint(4)+1 # input is True, generate
        else:
            print(' Not clear whether twirling; no twirling applied')
            self.twirl = None
            
        return self.twirl
        
    def _twirl(self):
        if not self.twirl: return # If twirl is passed None or False, do nothing
        i = self.twirl
        self.full_state = self.twirls[i-1] * self.full_state * self.twirls[i-1] # twirl matrices are Hermitian
    
    def _add_cOffset(self,constOffset,qubitDisplacement,add=True):
        '''
        Adds the cOffset parameter specified in mesolve_args to the random
        qubit offset specified by set_data_qubit_offsets
        
        By setting self.add_qubit_offsets = False this addition can be turned
        off and the qubit displacements are simply those specified by
        set_data_qubit_offsets()
        '''
        if add:
            return [sum(x) for x in zip(constOffset, qubitDisplacement)]
        else:
            return qubitDisplacement
    
    def _cOffset_cycle_correction(self,offset,cycle):
        '''
        Corrects the displacement offset to match the quadrant for
        the current cycle. Allows proper correlated errors due to
        
        '''
        if cycle%4==0:
            return [offset[0],offset[1],offset[2]]
        elif cycle%4==1:
            return [-offset[1],offset[0],offset[2]]
        elif cycle%4==2:
            return [-offset[0],-offset[1],offset[2]]
        elif cycle%4==3:
            return [offset[1],-offset[0],offset[2]]
        else:
            print('invalid cycle number '+str(cycle)+', unable to correct cOffset to match direction')
            return offset
    
    def run(self,time,steps,cycles=4):
        '''
        Run simulation for given time and steps
        '''
        self.last_run_quarter_cycle = []
        self.last_run_all = []
        cycle_time = time/cycles
        cycle_steps = steps/cycles
        
        if self.qubit_offsets_bool:
            cOffset = self.args['cOpts']['cOffset'] # store cOffset to reset after displacement
        #self.initial_states.append(self.full_state) # Optionally store initial states
        
        #TWIRLING
        if self.twirl:
            print(' Twirl '+str(self.twirl)+': '+self.twirls_str[self.twirl-1])
        self._twirl()
        

        
        for cycle in range(cycles):
            print(' Data qubit '+str(cycle+1))
        
            if self.qubit_offsets_bool: # Set qubit offsets as defined
                self.args['cOpts']['cOffset'] = self._add_cOffset(cOffset,self.qubit_offsets_list[cycle],self.add_qubit_offsets)
            # Correct displacement direction for cycle
            self.args['cOpts']['cOffset'] = self._cOffset_cycle_correction(self.args['cOpts']['cOffset'],cycle)
            
            # Print displacement if any is above 1 picometre
            if not np.allclose(self.args['cOpts']['cOffset'],0,atol=1e-12):
                print(' Data qubit displacement ' + str(self.args['cOpts']['cOffset']))
                # Displacement displayed in simulation reference frame, not lab frame
            
            probe_qubit = self.full_state.ptrace(0)
            data_qubit = self.full_state.ptrace(cycle+1)
            
            state = tensor(probe_qubit,data_qubit)
            state = self._run_quarter_cycle(state,cycle_time,cycle_steps)
            self.last_run_quarter_cycle.append(state)
            
            all_data_qubits = [self.full_state.ptrace(i+1) for i in range(4)]
            probe_qubit = state.ptrace(0)
            all_data_qubits[cycle] = state.ptrace(1)
            
            self.full_state = tensor(probe_qubit,all_data_qubits[0],all_data_qubits[1],all_data_qubits[2],all_data_qubits[3])
        
        if self.qubit_offsets_bool:
            self.args['cOpts']['cOffset'] = cOffset # restore cOffset in case it's important later
        
        #UNTWIRLING
        self._twirl()
        
        self.final_states.append(self.full_state)
        self._store_metadata(time,steps,cycles)
        
        return self.full_state