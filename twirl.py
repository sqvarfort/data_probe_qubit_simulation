from H_RWA import H_RWA
from Plotter import Plotter
from random import random
from DataQubitDisplacement import *
from Simulation import Simulation
from qutip import *
import yaml
import os

""" Simulation parameters """

info = os.path.join('twirl','twirl_config.yml')
config = os.path.join('twirl',"ham_conf.yml")
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
bitflips = lind_args.get('bit_errors')
initial_states = []
initial_states.append((pi/2.,0)) # probe qubit
for bit in bitflips: # data qubits
    if int(bit)==0:
        initial_states.append((0,0))
    elif int(bit)==1:
        initial_states.append((pi,0))
    else:
        print('Invalid bit_error specified')
print('Bitflips: '+str(bitflips))

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
displacements = lind_args.get('displacement')
displacements = [[float(x) for x in y] for y in displacements] # conversion from str to float
sim.set_data_qubit_offsets(displacements)

""" Set up folder/subfolder """
lind_args['folder'] = os.path.join(lind_args.get('folder'),lind_args.get('subfolder'))
if not os.path.exists(lind_args.get('folder')):
    os.mkdir(lind_args.get('folder'))
if not os.path.exists(os.path.join(lind_args.get('folder'),'data')):
    os.mkdir(os.path.join(lind_args.get('folder'),'data'))    


""" Run simulation """
for i in range(0,no_of_runs):
    print 'Starting loop ' + str(i+1)
    sim.choose_twirl(lind_args.get('twirl')) # Pass this False if you don't want twirling
    sim.run(time,steps)
    result_states = sim.last_run_all
    qsave(result_states, os.path.join(lind_args.get('folder'),'data','run'+str(i+1)))
    step_data = sim.last_run_quarter_cycle
    final_states.append(sim.last_run_all[-1])
    sim.reset_system_state()

qsave(final_states, os.path.join(lind_args.get('folder'),'data',"final_states"))
print('States saved')
for t in range(0,len(result_states),10):
    db.add_states(result_states[t].ptrace(0), kind='point')
    db.add_states(result_states[t].ptrace(1), kind='point')

db.save(os.path.join(lind_args.get('folder'),'bloch.pdf'))
#db.show()
print('Bloch for last run saved to '+os.path.join(lind_args.get('folder'),'bloch.pdf'))
Plotter(final_states, header,filetype='.pdf',display=False)
print('Plotter output to '+os.path.join(lind_args.get('folder'),''))