from H_RWA import H_RWA
from Plotter import Plotter

from DataQubitDisplacement import *
from Simulation import Simulation

from qutip import *
import yaml

Plotter.extract_phi()
""" Simulation parameters """

# here abrupt movement. tau is the time spend in interaction with each individual qubit

tau=7.7155254093495203e-05
info = 'info_config.yml'
config = "ham_conf.yml"
lind_args = dict()
config_args = dict()

if type(config) == str:
    with open(config) as cfile:
        config_args.update(yaml.load(cfile))
elif type(config) == dict:
    config_args.update(config)
else:
    print "config is wrong type"

if type(info) == str:
    with open(info) as cfile:
        lind_args.update(yaml.load(cfile))
elif type(config) == dict:
    lind_args.update(info)
else:
    print "info is wrong type"


header = lind_args.copy()
header.update(config_args)

""" Prepare Initial states """
intial_states = []
if lind_args.get('preparation_error'): # Check if preparation error should be applied
    if random.random() > 0.99: # Apply a preparation error at random
        initial_states.append((pi/2.,pi))
    else:
        intial_states.append((pi/2., 0))

if lind_args.get('odd'):
    initial_statesa.. = [(pi/2., pi), (0,0), (0,0), (0,0), (pi, 0)]




""" Prepare Hamiltonian """
ham=H_RWA(config)
hamiltonian = ham.getHfunc()

""" Prepare Bloch Sphere settings """
db=Bloch()
colors = ["g"]#,"r","g","#CC6600"]#,"r","g","#CC6600"]
db.point_color = "r"
db.point_marker = ['o']

""" Set up simulation """
mesolve_args = ham.getArgs()
iterations = 200
no_of_runs = 2
sim = Simulation(hamiltonian,initial_states,mesolve_args)
final_states = []


"""Adding Lindblad operators"""
if lind_args.get('dephasing'):
    sim.lind.dephasing(float(lind_args.get(dephasing_param)))
if lind_args.get('excitation'):
    sim.lind.excitation(float(lind_args.get(excitation_param)))
if lind_args.get('relaxation'):
    sim.lind.relaxation(float(lind_args.get(relaxation_param)))



""" Run simulation """
for i in range(1,no_of_runs):
    print 'Starting loop ' + str(i)
    sim.run(tau,iterations,1)
    result_states = sim.last_run_all
    step_data = sim.last_run_quarter_cycle
    final_states.append(sim.last_run_all[-1])
    sim.reset_system_state()



for t in range(0,len(result_states),10):
    db.add_states(result_states[t].ptrace(0), kind='point')
    db.add_states(result_states[t].ptrace(1), kind='point')
db.show()

Plotter(final_states, header)
