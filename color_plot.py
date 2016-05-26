from H_RWA import H_RWA
from Plotter import Plotter
from random import random
from DataQubitDisplacement import *
from Simulation import Simulation
from qutip import *
import yaml

""" Prepare Bloch Sphere settings """
db=Bloch()
colors = ["b","g","r","k"]#,"r","g","#CC6600"]
db.point_color = colors
db.point_marker = ['o']


result_states = qload("anim/notodd")

for i in range(4):
	db.add_points(np.transpose(np.array([Plotter.bloch_vector(result_states[t].ptrace(0)) for t in range(int(i*config_args['steps']/4.),int((i+1)*config_args['steps']/4.),10)])))

db.show()