#Plotting functions

from numpy import *
import matplotlib.pyplot as pltB
from qutip import *
from pylab import *
from scipy import constants as cp
import os



class Plotter(object):
    """ Class that calculates statistics and plots phase on histograms"""

    def __init__(self, states, info):
        self.states= states
        self.info = info
        # Extract list of probe states
        probe_states = [states[t].ptrace(0) for t in range(0, len(states))]
        list_of_phi = [self.extract_phi(self.bloch_vector(probe_states[t])) for t in range(0,len(probe_states))]
        self.Histogram(list_of_phi)
        self.save_to_file(list_of_phi)

    def bloch_vector(self, state): # Calculate the bloch vector
        return [expect(sigmax(), state),
                   expect(sigmay(), state),
                   expect(sigmaz(), state)]

    def extract_phi(self, vec):
        if vec[0] < 0:
            return arctan(vec[1]/vec[0])
        else:
            return arctan(vec[1]/vec[0]) + pi


    def Histogram(self, list_of_phi):


        plt.hist(list_of_phi)
        plt.show()
        ax = plt.axes()
        plt.savefig('Histogram.pdf')

    def save_to_file(self, list_of_phi, info):
        os.chdir(path) # Set path where rates are saved to current directory
        f = open('phi_raw_data.txt', 'w')
        res.write('Values of phi \n')
        res.write('Lindblad operators:' + info)
        res.write('Number of runs:' + len(list_of_phi))
        # Write all values to the file
        for i in len(list_of_phi):
            res.write(list_of_phi[i] + "\n")
        # Calculate mean and std and write them to file
        mean_phi = mean(list_of_phi)
        std_phi = std(list_of_phi)
        res.write('The mean phase is:' + mean_phi + 'which is' + abs(mean_phi - pi/2.)/(pi/2.) + "% of pi/2" + "\n")
        res.write("\n" + 'The standard deviation is:' + std_phi + "which is" + abs(std_phi - pi/2.)/(pi/2.) + "% of pi/2" + "\n")





#def Histogram()
