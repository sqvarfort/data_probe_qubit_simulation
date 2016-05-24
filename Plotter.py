#Plotting functions

from numpy import *
import matplotlib.pyplot as pltB
from qutip import *
from pylab import *
from scipy import constants as cp
import os
from matplotlib import rc
import json



class Plotter(object):
    """ Class that calculates statistics and plots phase on histograms
        Input:      - States: An array of states to be plotted an processed.
                    - Info: Information about current simulation, goes in the header of the files
        Output:     - histogram (pdf)
                    - textfile of all phi values
                    - textfile of all probe-data full states
                    - textfile of all probe-qubit states
                    - textfile of all data-qubit states
    """

    def __init__(self, states, info):
        self.states = states
        self.info = info

        # Save the states in the array to a file
        #self.states_to_file(states, info)

        # Extract list of probe states
        probe_states = [states[t].ptrace(0) for t in range(0, len(states))]

        # Calculate phi for every state
        list_of_phi = [self.extract_phi(self.bloch_vector(probe_states[t])) for t in range(0,len(probe_states))]

        # Plot histogram and write textfile
        self.Histogram(list_of_phi)
        self.phi_to_file(list_of_phi, info)


    def bloch_vector(self, state): # Calculate the bloch vector
        return [expect(sigmax(), state),
                   expect(sigmay(), state),
                   expect(sigmaz(), state)]

    def extract_phi(self, vec): # Extract values for phi. If-statements needed for stupid non-invertible arctan
        if vec[0] >= 0.:
            if vec[1] < 0.:
                return -arctan(vec[1]/vec[0])
            else:
                return 2.*pi - arctan(vec[1]/vec[0])
        else:
            return pi - arctan(vec[1]/vec[0])

    def Histogram(self, list_of_phi):
        # Use LaTeX rendering
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        path = os.path.abspath(__file__) # Set path where rates are saved to current directory
        fig = plt.figure()
        ax = fig.add_subplot(111,title='Phase Histogram')
        ax.set_xlabel('$\phi$')
        plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
           ['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']) #Trig ticks
        ax.set_ylabel('N')
        #plt.plot(list_of_phi)
        plt.hist(list_of_phi, bins = 20, label = 'No. of runs:' + str(len(list_of_phi)))
        ax.legend()
        ax.legend(loc='upper right')
        plt.show()
        plt.savefig('Histogram.pdf')

    def phi_to_file(self, list_of_phi, info):
        path = os.path.abspath(__file__)
        f = open('phi_data.txt', 'w+')
        with open('phi_data.txt', 'w+') as fp:
            json.dump(info, fp)

        # Write all values to the file
        for item in list_of_phi:
            f.write("%s\n" % item)

        # Calculate mean and std and write them to file
        mean_phi = self.calc_mean(list_of_phi)
        std_phi = self.calc_std(list_of_phi)
        f.write('The mean phase is: ' + str(mean_phi) + ' which is ' + str(abs(mean_phi - pi/2.)/(pi/2.)) + "% of pi/2" + "\n")
        f.write("\n" + 'The standard deviation is: ' + str(std_phi) + " which is " + str(abs(std_phi - pi/2.)/(pi/2.)) + "% of pi/2" + "\n")

    def calc_mean(self, list_of_phi):
        return mean(list_of_phi)

    def calc_std(self, list_of_phi):
        return std(list_of_phi)

    def states_to_file(self, states, info):
        probe_states = [states[t].ptrace(0) for t in range(0, len(states))]
        data_states = [states[t].ptrace(1) for t in range(0, len(states))]
        path = os.path.abspath(__file__)
        fprobe_states = open('probe_states.txt', 'w+')
        fdata_states = open('data_states.txt', 'w+')
        fprobe_data = open('full_states.txt', 'w+')

        fprobe_data.write('Full probe and data states \n')
        fprobe_states.write('Traced out probe states \n')
        fprobe_states.write('Traced out probe states \n')

        fprobe_data.write(info + '\n') # Write info into the
        fprobe_states.write(info + '\n') # Write info into the
        fdata_states.write( info + '\n') # Write info into the

        states_output = np.vstack((states))
        probe_states_output = np.vstack((probe_states))
        data_states_output = np.vstack((data_states))

        file_data_store('full_states.txt', states_output, numtype="complex")
        file_data_store('probe_states.txt', probe_states_output, numtype = "complex")
        file_data_store('data_states.txt', data_states_output, numtype = "complex")


    def Probe_measurement_outcome(self, probe_states):
        """ This method calculates the probability of measuring the probe qubit in a particular state after one round of the measurement. We expect the probe qubit to either end up in the |+> state for an even outcome, or in the |-> state for the odd outcome.
        """
        fprobes = open('probe_measurements.txt', 'w+')
        fprobes.write('Probe Measurement Outcomes \n')
        fprobes.write(info + '\n')
        fprobes.write('sigmax exp.' + '\t' + 'Prob. of |+>'+ '\t'+ 'Prob. of |->' + '\n')
        measurements = [expect(sigmax(), state) for state in probe_states]

        # Given the expectation value for a sigmax measurement, we solve for the probailities of measuring |+> and |->

        plus_prob = [((1/2.)*(E + sqrt(2-E**2))) for E in measurements]
        min_prob = [(1/2.)*(-E + sqrt(2-E**2)) for E in measurements]
        avg_exp = mean(measurements)
        plus_mean = mean(plus_prob)
        min_mean = mean(min_prob)

        for Exp in measurements:
            fprobes.write("%s\t %s \t %s \n" % entry, (1/2.)*(Exp + sqrt(2-Exp**2)), (1/2.)*(-Exp + sqrt(2-Exp**2)))

        fprobes.write('The average expectation value is: ' + str(avg_exp) + '\n')
        fprobes.write('The average probability of measuring |+> is: ' + str(plus_mean) + '\n')
        fprobes.write('The average probability of measuring |-> is: ' + str(min_mean) + '\n')











#def Histogram()
