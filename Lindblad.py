from math import sqrt
import numpy as np
from qutip import *

class Lindblad():
    def __init__(self,qubits=1):
        self.lindblads = []
        self.qubits = qubits
        
    def clear(self):
        self.lindblads = []
        
    def _prepare_tensor(self,lindblad,qubit):
        qubit += 1 # Account for addressing 1st qubit with 0
        qubits_before = qubit-1
        qubits_after = self.qubits-qubit
        
        return tensor(qeye(qubits_before+1),lindblad,qeye(qubits_after+1))
        
    def relaxation(self,rate,qubit=0):
        mat = np.matrix([[0,1],[0,0]])
        single_lindblad = Qobj(sqrt(rate)*mat)
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
    def excitation(self,rate,qubit=0):
        mat = np.matrix([[0,0],[1,0]])
        single_lindblad = Qobj(sqrt(rate)*mat)
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
    def dephasing(self,rate,qubit=0):
        mat = np.matrix([[1,0],[0,-1]])
        single_lindblad = Qobj(sqrt(0.5*rate)*mat) # Factor of 0.5 makes rate match exponential decay
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
        
    #------ Below this line things get less physically meaningful ------#
        
    def _depolarising(self,rate,qubit=0): # does not give expected behaviour
        I = np.matrix([[1,0],[0,1]])
        x = np.matrix([[0,1],[1,0]])
        y = np.matrix([[0,-1j],[1j,0]])
        z = np.matrix([[1,0],[0,-1]])
        ind_rate = sqrt(rate/4.0) # sqrt(p'/3) where p'/p*(4/3)
        single_lindblad = Qobj( sqrt(1-rate)*I +  ind_rate*x + ind_rate*y + ind_rate*z )
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
    def sigmax(self,rate,qubit=0): # Not sure this is so physically meaningful
        mat = np.matrix([[0,1],[1,0]])
        single_lindblad = Qobj(sqrt(rate)*mat)
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
    def sigmay(self,rate,qubit=0): # Not sure this is so physically meaningful
        mat = np.matrix([[0,-1j],[1j,0]])
        single_lindblad = Qobj(sqrt(rate)*mat)
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
    def sigmaz(self,rate,qubit=0): # Basically same as dephasing
        mat = np.matrix([[1,0],[0,-1]])
        single_lindblad = Qobj(sqrt(rate)*mat)
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad
        
    def identity(self,rate,qubit=0): # Does nothing with rate p
        mat = np.matrix([[1,0],[0,1]])
        single_lindblad = Qobj(sqrt(rate)*mat)
        overall_lindblad = self._prepare_tensor(single_lindblad,qubit)
        self.lindblads.append(overall_lindblad)
        return overall_lindblad