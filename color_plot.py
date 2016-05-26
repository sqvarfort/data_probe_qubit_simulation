from H_RWA import H_RWA
from Plotter import Plotter
from random import random
from DataQubitDisplacement import *
from Simulation import Simulation
from qutip import *
import yaml
import numpy as np

config = "ham_conf.yml"
config_args = {}

if type(config) == str:
    with open(config) as cfile:
        config_args.update(yaml.load(cfile))
elif type(config) == dict:
    config_args.update(config)
else:
    print "config is wrong type"


""" Prepare Bloch Sphere settings """
db=Bloch()
colors = ["b","g","r","k"]#,"r","g","#CC6600"]
db.point_color = colors
db.point_marker = ['o']


result_states = qload("anim/odd")

for i in range(4):
	db.add_points(np.transpose(np.array([Plotter.bloch_vector(result_states[t].ptrace(0)) for t in range(int(i*config_args['steps']/4.),int((i+1)*config_args['steps']/4.),10)])))

db.show()
