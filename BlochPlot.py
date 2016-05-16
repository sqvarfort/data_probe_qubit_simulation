import numpy as np
import scipy as sp
import qutip as qt

def bloch_plot(dm_array,kind='point'):
	'''
	Given an array of single-qubit density matrices,
	plots density matrices on Bloch sphere
	and displays using matplotlib
	
	Second arg optionally specifies whether 'point' or 'vector'
	'''
	b = qt.Bloch()
	
	if type(dm_array) == list: # Input is array of density matrices to be plotted
		qobj_array = []
		for dm in dm_array:
			qobj_array.append( qt.Qobj(dm) )
		b.add_states(qobj_array,kind)
		b.show()
		
	else: # Assume input is one density matrix
		qobj = qt.Qobj(dm_array)
		b.add_states(qobj)
		b.show()
		
		
# TODO: Make points plotted at same time the same colour