from H_RWA import H_RWA
from Plotter import Plotter
from random import random
from DataQubitDisplacement import *
from Simulation import Simulation
from qutip import *
import yaml
import os
import multiprocessing as mp
import time


""" Simulation parameters """

# here abrupt movement. tau is the time spend in interaction with each individual qubit

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
initial_states = []
if lind_args.get('preparation_error'): # Check if preparation error should be applied
    if random.random() > 0.99: # Apply a preparation error at random
        initial_states.append((pi/2.,pi))
    else:
        initial_states.append((pi/2., 0))
else:
    initial_states.append((pi/2.,0))
if lind_args.get('odd'):
    for i in range(0,3):
        initial_states.append((0,0))
    initial_states.append((pi, 0))
else:
    for i in range(0,4):
        initial_states.append((0,0))


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
steps = config_args['steps']
time = config_args['time']
no_of_runs = lind_args['runs']
sim = Simulation(hamiltonian,initial_states,mesolve_args)
final_states = []

"""Adding Lindblad operators"""
if lind_args.get('dephasing'):
    sim.lind.dephasing(lind_args.get('dephasing_param'))
if lind_args.get('excitation'):
    sim.lind.excitation(lind_args.get('excitation_param'))
if lind_args.get('relaxation'):
    sim.lind.relaxation(lind_args.get('relaxation_param'))

""" Set up and generate qubit displacements """
if lind_args.get('qubit_displacement_error'):
    disp_radius = float(lind_args.get('qubit_displacement_radius'))
    disp_halfheight = float(lind_args.get('qubit_displacement_halfheight'))

    sim.generate_data_qubit_offsets(disp_radius,disp_halfheight)


def log_result(result):
    # This is called whenever simulation(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    final_states.insert(result[0], result[1])

def do_simulation(i):
    print 'Starting loop ' + str(i)
    sim.regenerate_data_qubit_offsets() # Necessary for each run to have different displacement errors
	#sim.args['J'] = jlist[i]
    sim.run(time/4,steps,1)	#adjust for full circle!
    result_states = sim.last_run_all
    qsave(result_states, os.path.join(lind_args['folder'],"iteration"+str(i)))
    #step_data = sim.last_run_quarter_cycle
    sim.reset_system_state()
    # return tuple of iteration and result, because jobs are started async = return whenever finished
    return (i, sim.last_run_all[-1]) 


if __name__ == '__main__': # windows needs this for multiprocessing ... whatever
	""" Run simulation """
	pool = mp.Pool(processes=8)# limits no of parallel processes to 3 = 75% CPU usage for me. more slows laptop pc down
	# windows pc: 1 process=120s each. 2 processes= 130s each (nearly twice as fast). 
	# 3processes= 130s each (3 times as fast) 4 processes= 138s, 6 processes=210s each
	# still speed up for 10: 337s each. 
	mp.freeze_support()
	for i in range(0,no_of_runs):
		# we could parse args in addition to i and change significant parameters during interations!
		pool.apply_async(do_simulation, args = (i, ), callback = log_result) 
		
	pool.close()
	pool.join()
	#print final_states
	qsave(final_states, os.path.join(lind_args.get('folder'),"final_states"))

"""
for t in range(0,len(result_states),10):
    db.add_states(result_states[t].ptrace(0), kind='point')
    db.add_states(result_states[t].ptrace(1), kind='point')
db.show()

Plotter(final_states, header)
"""
