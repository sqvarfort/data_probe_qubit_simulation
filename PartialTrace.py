import numpy as np
import scipy as sp
import qutip as qt

def partial_trace(dm,qubit_index=0):
	'''
	For 4x4 density matrix dm, returns the partial trace
	leaving the qubit indexed by qubit_index
	0: First qubit in Kronecker product
	1: Second qubit in Kronecker product
	'''
	
	Qobj = qt.Qobj(dm,[[2,2], [2,2]]) # Initialise a Qobj representing dm of 2 qubits
	Qobj = Qobj.ptrace(qubit_index)
	output_dm = Qobj.data.todense() # Retrieve Qobj data as np.matrix
	
	return output_dm
	
	
def decompose(dm):
	'''
	For 4x4 density matrix dm, returns the partial traces for
	both qubits as 2x2 matrices in an array
	'''
	pt_probe = partial_trace(dm,0)
	pt_data = partial_trace(dm,1)
	
	return [pt_probe,pt_data]
	
	
	
# Some example code that traces out 4x4 matrix dm to retrieve original matrices
# dm_a = np.array([[1,0],[0,0]])
# dm_b = np.array([[0,0],[0,1]])

# dm = np.kron(dm_a,dm_b)

# print decompose(dm)