from H_RWA import H_RWA
from Plotter import Plotter
from random import random
from DataQubitDisplacement import *
from Simulation import Simulation
from qutip import *
import yaml
import os, shutil
import multiprocessing as mp
import time

""" Simulation parameters """

info = os.path.join('displacement','displacement_config.yml')
config = os.path.join('displacement',"ham_conf.yml")
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


#header = lind_args.copy()
#header.update(config_args)

""" Prepare Initial states """
bitflips = lind_args.get('bit_errors')

def calculate_initial_states_from_bitflips(bitflips):
    initial_states = []
    initial_states.append((pi/2.,0)) # probe qubit
    for bit in bitflips: # data qubits
        if int(bit)==0:
            initial_states.append((0,0))
        elif int(bit)==1:
            initial_states.append((pi,0))
        else:
            print('Invalid bit_error specified')
            initial_states.append((0,0))
    print('Bitflips: '+str(bitflips))
    return initial_states

#initial_states = calculate_initial_states_from_bitflips(bitflips)
    
""" Prepare Hamiltonian """
ham=H_RWA(config)
hamiltonian = ham.getHfunc()

""" Prepare Bloch Sphere settings """
def prepare_bloch():
    db=Bloch()
    colors = ["g"]#,"r","g","#CC6600"]#,"r","g","#CC6600"]
    db.point_color = "r"
    db.point_marker = ['o']
    return db

def plot_bloch(db,result_states):
    for t in range(0,len(result_states),10):
        db.add_states(result_states[t].ptrace(0), kind='point')
        db.add_states(result_states[t].ptrace(1), kind='point')    
    
""" Saving twirls applied """
displacement_applied = []
def save_displacement(displacement_applied):
    f = open(os.path.join(lind_args.get('folder'),lind_args.get('subfolder'),'data','displacements.txt'),'w')
    f.write( yaml.dump(displacement_applied))
    f.close()

""" Set up simulation """
mesolve_args = ham.getArgs()
steps = config_args['steps']
time = config_args['time']
no_of_runs = lind_args['runs']
# sim = Simulation(hamiltonian,initial_states,mesolve_args) # Need new simulation for each run
final_states = []


# """Adding Lindblad operators"""
# if lind_args.get('dephasing'):
    # sim.lind.dephasing(lind_args.get('dephasing_param'))
# if lind_args.get('excitation'):
    # sim.lind.excitation(lind_args.get('excitation_param'))
# if lind_args.get('relaxation'):
    # sim.lind.relaxation(lind_args.get('relaxation_param'))

""" Set up and generate qubit displacements """
#displacements = lind_args.get('displacement')
#displacements = [[float(x) for x in y] for y in displacements] # conversion from str to float
#sim.set_data_qubit_offsets(displacements)

""" Set up and generate qubit displacements """
disp_radius = float(lind_args.get('qubit_displacement_radius'))
disp_halfheight = float(lind_args.get('qubit_displacement_halfheight'))




""" Set up folder/subfolder and save header """
#lind_args['folder'] = os.path.join(lind_args.get('folder'),lind_args.get('subfolder'))
if not os.path.exists(lind_args.get('folder')):
    os.mkdir(lind_args.get('folder'))
if not os.path.exists(os.path.join(lind_args.get('folder'),lind_args.get('subfolder'))):
    os.mkdir(os.path.join(lind_args.get('folder'),lind_args.get('subfolder')))
if not os.path.exists(os.path.join(lind_args.get('folder'),lind_args.get('subfolder'),'data')):
    os.mkdir(os.path.join(lind_args.get('folder'),lind_args.get('subfolder'),'data')) 
    
header = lind_args.copy()
header.update(config_args)
def save_header(header):
    f = open(os.path.join(lind_args.get('folder'),lind_args.get('subfolder'),'data','header.txt'),'w')
    f.write( yaml.dump(header))
    f.close()
save_header(header)


""" Copy input .yml files to subfolder """
shutil.copy2(info,os.path.join(lind_args.get('folder'),lind_args.get('subfolder')))
shutil.copy2(config,os.path.join(lind_args.get('folder'),lind_args.get('subfolder')))

    
def log_result(result):
    # This is called whenever simulation(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    final_states.insert(result[0], result[1])
    displacement_applied.insert(result[0], result[2])

def do_simulation(i):
    print 'Starting loop ' + str(i)
    
    if lind_args.get('specify_errors'): initial_states = calculate_initial_states_from_bitflips(bitflips[i])
    else: initial_states = calculate_initial_states_from_bitflips(bitflips)
    
    sim = Simulation(hamiltonian,initial_states,mesolve_args)
    print '1'
    sim.generate_data_qubit_offsets(disp_radius,disp_halfheight)
    print '2'
    sim.loop = i
    #sim.set_data_qubit_offsets(displacements)
    #sim.choose_twirl(lind_args.get('twirl'))
    sim.run(time,steps)    #adjust for full circle!
    print '3'
    result_states = sim.last_run_all
    qsave(result_states, os.path.join(lind_args.get('folder'),lind_args.get('subfolder'),'data','run'+str(i+1)))
    #step_data = sim.last_run_quarter_cycle
    sim.reset_system_state()
    # return tuple of iteration and result, because jobs are started async = return whenever finished
    return (i, sim.last_run_all[-1], sim.qubit_offsets_list) 


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
    qsave(final_states, os.path.join(lind_args.get('folder'),lind_args.get('subfolder'),'data',"final_states"))
    save_displacement(displacement_applied)
    
    

# """ Run simulation """
# for i in range(0,no_of_runs):
    # print 'Starting loop ' + str(i+1)
    # sim.choose_twirl(lind_args.get('twirl')) # Pass this False if you don't want twirling
    # sim.run(time,steps)
    # result_states = sim.last_run_all
    # qsave(result_states, os.path.join(lind_args.get('folder'),'data','run'+str(i+1)))
    # step_data = sim.last_run_quarter_cycle
    # final_states.append(sim.last_run_all[-1])
    # sim.reset_system_state()

# qsave(final_states, os.path.join(lind_args.get('folder'),'data',"final_states"))
# print('States saved')
# for t in range(0,len(result_states),10):
    # db.add_states(result_states[t].ptrace(0), kind='point')
    # db.add_states(result_states[t].ptrace(1), kind='point')

# db.save(os.path.join(lind_args.get('folder'),'bloch.pdf'))
# db.show()
# print('Bloch for last run saved to '+os.path.join(lind_args.get('folder'),'bloch.pdf'))
# Plotter(final_states, header,filetype='.pdf',display=False)
# print('Plotter output to '+os.path.join(lind_args.get('folder'),''))