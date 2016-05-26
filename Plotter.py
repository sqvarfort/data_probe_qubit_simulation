#Plotting functions

from numpy import *
import matplotlib.pyplot as plt
from qutip import *
from pylab import *
from scipy import constants as cp
import os
from matplotlib import rc
import yaml
import datetime
import time




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

    def __init__(self, states, info,filetype='.png',display=True):
        self.states = states
        self.info = info
        self.filetype = filetype
        self.display = display

        # Save the states in the array to a file
        #self.states_to_file(states, info)

        # Extract list of probe states
        probe_states = [states[t].ptrace(0) for t in range(0, len(states))]

        # Calculate phi for every state
        list_of_phi = [self.extract_phi(self.bloch_vector(probe_states[t])) for t in range(0,len(probe_states))]

        # Plot histogram and write textfile

        self.Histogram(list_of_phi, self.info)
        self.phi_to_file(list_of_phi, self.info)
        self.Probe_measurement_outcome(probe_states, self.info)
        self.states_to_file(states, self.info)



    @staticmethod
    def bloch_vector(state): # Calculate the bloch vector
        return [expect(sigmax(), state),
                   expect(sigmay(), state),
                   expect(sigmaz(), state)]
    @staticmethod
    def extract_phi(vec): # Extract values for phi. If-statements needed for stupid non-invertible arctan
        if vec[0] >= 0.:
            if vec[1] <= 0.:
                return -arctan(vec[1]/vec[0])
            else:
                return 2.*pi - arctan(vec[1]/vec[0])
        else:
            return pi - arctan(vec[1]/vec[0])

    def Histogram(self, list_of_phi, info):
        # Use LaTeX rendering
        folder_name = info.get('folder')
        if 'subfolder' in info:
            folder_name = os.path.join(info.get('folder'),info.get('subfolder'))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel('$\phi$')
        plt.ylabel('N')
        plt.title('Phase Histogram')
        plt.grid(True)

        # Move data-points that are close to zero to other end
        for item in list_of_phi:
            if item < pi/2.:
                item = item + 2.*pi

        n, bins, patches = plt.hist(list_of_phi, 100, facecolor='green', alpha=0.75, range = (0,2*pi), label = 'No. of runs: ' + str(len(list_of_phi)))


        #Ticks go from pi/2 to pi/2 including the origin
        plt.xticks([pi/2, pi, 3*pi/2, 2*pi, 5*pi/2, 3*pi], [r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$0 (2\pi)$', r'$\frac{\pi}{2}$'])


        # If there is an error in the data qubits,
        """
        if info_config.get('odd'):
            plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi], ['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
        else:
            for item in list_of_phi:
                if item < pi:
                    item = item + 2.*pi
            plt.xticks([pi, 3*pi/2, 2*pi, 5*pi/2, 3*pi], [r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'], '$0$', r'$\frac{\pi}{2}$')
        """

        # Render plot
        plt.legend()
        plt.legend(loc='upper right')
        fig = plt.gcf()
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')
        plt.savefig(os.path.join(folder_name,'Histogram') + st + self.filetype,  transparent=False, dpi=1000)
        if self.display:
            plt.show()



    def phi_to_file(self, list_of_phi, info):
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')
        folder_name = info.get('folder')
        if 'subfolder' in info:
            folder_name = os.path.join(info.get('folder'),info.get('subfolder'))
        f = open(os.path.join(folder_name,'phi_data') + st + '.txt', 'w+')
        f.write('Phase of the final state' '\n \n')
        f.write( yaml.dump(info, default_flow_style=False) + '\n \n')
        # Write all values to the file
        for item in list_of_phi:
            f.write("%s\n" % item)

        # Calculate mean and std and write them to file
        mean_phi = self.calc_mean(list_of_phi)
        std_phi = self.calc_std(list_of_phi)
        f.write('The mean phase is: ' + str(mean_phi))
        f.write( 'The mean error is ' + str(100.*abs(mean_phi - 2.*pi)/(2.*pi)) + "% of 2pi, and " + str(100.*abs(mean_phi - pi)/pi) + '% of pi' "\n")
        f.write("\n" + 'The standard deviation is: ' + str(std_phi) + " which is " + str(100.*abs(std_phi - 2.*pi)/(2.*pi)) + "% of 2pi and " + str(100.*abs(std_phi - pi)/pi) + "% of pi." "\n")


    def calc_mean(self, list_of_phi):
        return mean(list_of_phi)

    def calc_std(self, list_of_phi):
        return std(list_of_phi)

    def states_to_file(self, states, info):
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')
        folder_name = info.get('folder')
        if 'subfolder' in info:
            folder_name = os.path.join(info.get('folder'),info.get('subfolder'))
        probe_states = [states[t].ptrace(0) for t in range(0, len(states))]
        data_states = [states[t].ptrace(1) for t in range(0, len(states))]
        path = os.path.abspath(__file__)
        fprobe_states = open(os.path.join(folder_name,'probe_states') + st + '.txt', 'w+')
        fdata_states = open(os.path.join(folder_name,'data_states') + st + '.txt', 'w+')
        fprobe_data = open(os.path.join(folder_name,'full_states') + st + '.txt', 'w+')

        fprobe_data.write('Full probe and data states \n \n')
        fprobe_states.write('Traced out probe states \n \n')
        fdata_states.write('Traced out data states \n \n')

        fprobe_data.write( yaml.dump(info, default_flow_style=False) + '\n \n')
        fprobe_states.write( yaml.dump(info, default_flow_style=False) +'\n \n' )
        fdata_states.write( yaml.dump(info, default_flow_style=False) + '\n \n')

        states_output = np.vstack((states))
        probe_states_output = np.vstack((probe_states))
        data_states_output = np.vstack((data_states))

        for item in states_output:
            fprobe_data.write("%s\n" % item)

        for item in probe_states_output:
            fprobe_states.write("%s\n" % item)

        for item in data_states_output:
            fdata_states.write("%s\n" % item)


    def Probe_measurement_outcome(self, probe_states, info):
        """ This method calculates the probability of measuring the probe qubit in a particular state after one round of the measurement. We expect the probe qubit to either end up in the |+> state for an even outcome, or in the |-> state for the odd outcome.
        """
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')
        folder_name = info.get('folder')
        if 'subfolder' in info:
            folder_name = os.path.join(info.get('folder'),info.get('subfolder'))
        fprobes = open(os.path.join(folder_name,'probe_measurements') + st + '.txt', 'w+')
        fprobes.write('Probe Measurement Outcomes \n')

        fprobes.write( yaml.dump(info, default_flow_style=False) + '\n \n')

        fprobes.write('sigmax exp.' + '\t' + 'Prob. of |+>'+ '\t'+ 'Prob. of |->' + '\n')
        measurements = [expect(sigmax(), state) for state in probe_states]

        # Given the expectation value for a sigmax measurement, we solve for the probailities of measuring |+> and |->

        plus_prob = [((1/2.)*(E + sqrt(2-E**2))) for E in measurements]
        min_prob = [(1. -p) for p in plus_prob]
        avg_exp = mean(measurements)
        plus_mean = mean(plus_prob)
        min_mean = mean(min_prob)

        for Exp in measurements:
            fprobes.write("%s\t %s \t %s \n" % (Exp, ((Exp+1.)/2.), (1. - ((Exp+1.)/2.))))

        fprobes.write('The average expectation value is: ' + str(avg_exp) + '\n')
        fprobes.write('The average probability of measuring |+> is: ' + str(plus_mean) + '\n')
        fprobes.write('The average probability of measuring |-> is: ' + str(min_mean) + '\n')
