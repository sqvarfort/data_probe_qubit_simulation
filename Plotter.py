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
        self.states = states
        self.info = info
        # Extract list of probe states
        probe_states = [states[t].ptrace(0) for t in range(0, len(states))]
        list_of_phi = [self.extract_phi(self.bloch_vector(probe_states[t])) for t in range(0,len(probe_states))]
        self.Histogram(list_of_phi)
        self.save_to_file(list_of_phi, info)

    def bloch_vector(self, state): # Calculate the bloch vector
        return [expect(sigmax(), state),
                   expect(sigmay(), state),
                   expect(sigmaz(), state)]

    def extract_phi(self, vec):
        if vec[0] > 0.:
            if vec[1] < 0.:
                return -arctan(vec[1]/vec[0])
            else:
                return 2.*pi - arctan(vec[1]/vec[0])
        else:
            return pi - arctan(vec[1]/vec[0])

    def Histogram(self, list_of_phi):
        path = os.path.abspath(__file__) # Set path where rates are saved to current directory
        fig = plt.figure()
        ax = fig.add_subplot(111,title='Phase Histogram')
        ax.set_xlabel('phi')
        ax.set_ylabel('N')
        #plt.plot(list_of_phi)
        plt.hist(list_of_phi)
        plt.show()
        plt.savefig('Histogram.pdf', path)

    def save_to_file(self, list_of_phi, info):
        f = open('phi_raw_data.txt', 'w')
        res.write('Values of phi \n')
        res.write('Lindblad operators:' + info)
        res.write('Number of runs:' + len(list_of_phi))
        # Write all values to the file
        for i in range(0, len(list_of_phi)):
            res.write(list_of_phi[i] + "\n")

        # Calculate mean and std and write them to file
        mean_phi = self.calc_mean(list_of_phi)
        std_phi = self.calc_std(list_of_phi)
        res.write('The mean phase is:' + mean_phi + 'which is' + abs(mean_phi - pi/2.)/(pi/2.) + "% of pi/2" + "\n")
        res.write("\n" + 'The standard deviation is:' + std_phi + "which is" + abs(std_phi - pi/2.)/(pi/2.) + "% of pi/2" + "\n")

    def calc_mean(self, list_of_phi):
        return mean(list_of_phi)

    def calc_std(self, list_of_phi):
        return std(list_of_phi)



#def Histogram()
