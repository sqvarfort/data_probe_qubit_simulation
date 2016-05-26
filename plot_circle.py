from scipy import constants as cp
from pylab import *
from qutip import *
import yaml
from H_RWA import H_RWA

config = "ham_conf.yml"
cargs = dict()

if type(config) == str:
    with open(config) as cfile:
        cargs.update(yaml.load(cfile))
elif type(config) == dict:
    cargs.update(config)
else:
    print "config is wrong type"


tlist=linspace(0, cargs['time'], cargs['steps'])
circ=H_RWA.circ_motion(tlist, cargs['cOpts'])

plot(circ[0], circ[1])

#plot(tlist, circ[2])

show()

