#Phase plotting

from numpy import *
from matplotlib import pyplot as plt
import time
import datetime
from matplotlib import rcParams
from qutip import *
import yaml
from Plotter import Plotter
from scipy.optimize import curve_fit
from scipy.stats import exponpow
from matplotlib.patches import BoxStyle

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



#states = qload('dephasing/phys_params_circ/final_states')
#Plotter(states, header)


circ_array = loadtxt("dephasing/phys_params_circ/circ_phi_data.txt", comments="#", delimiter=",", unpack=False)

abrupt_array = loadtxt("dephasing/phys_params_abrupt/abrupt_phi_data.txt", comments="#", delimiter=",", unpack=False)

abrupt_probs = []
circ_probs = []
abrupt_dephasing = []
circ_dephasing = []

for i in range(0,len(abrupt_array)):
    abrupt_dephasing.append(abrupt_array[i][0])
    abrupt_probs.append(abrupt_array[i][1])



for i in range(0,len(circ_array)):
    circ_dephasing.append(circ_array[i][0])
    circ_probs.append(circ_array[i][1])

fig = plt.figure()




# Order data-points
"""
circ_dephasing = circ_dephasing[circ_dephasing.argsort()]
circ_probs = circ_probs[circ_probs.argsort()]
abrupt_dephasing = abrupt_dephasing[abrupt_dephasing.argsort()]
abrupt_probs = abrupt_probs[abrupt_probs.argsort()]
"""


# Make x-data logarithmic
circ_dephasing = log10(circ_dephasing)
abrupt_dephasing = log10(abrupt_dephasing)







# Fit lines
def func(x, a, b, c):
    return a * exp(-b * x) + c

# Define doubly exponential function

def func1(x, a):
    return a*exp(-a*x)

#a*x + b*x**2 + c*x**3 + d
#(a/b)*(x/b)**(a-1)*exp(-(x/b)**a)
#1./((1/2.)*exp((x-a)/b) + c)
#(exp(1-exp(a*x**b))*exp(a*x**b)*a*b*x**(b-1))


plt.plot(circ_dephasing, circ_probs, '-o', label = 'Circular orbit', alpha = 0.5, markersize=35, linewidth=2)
plt.plot(abrupt_dephasing, abrupt_probs, '-o', label = 'Abrupt orbit', alpha = 0.5, markersize=35, linewidth=2)

#plt.plot(linspace(-1.,10.0,2000), func(linspace(-1.0,10.0,2000), popt1[0], popt1[1], popt1[2]), 'r-', linewidth=2)

#plt.plot(linspace(0,10.0,2000),0.5 + func1(linspace(0,10.0,2000),*popt2), 'r-')


for i in range(4, 9):
    plt.plot((circ_dephasing[i], abrupt_dephasing[i]), (circ_probs[i], abrupt_probs[i]), 'k-', linewidth = 2)

#fig_width_pt = 276  # Get this from LaTeX using \showthe\linewidth
#inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
#fig_width = fig_width_pt*inches_per_pt  # width in inches
#fig_height = fig_width*golden_mean      # height in inches
#fig_size =  [fig_width,fig_height]#in inch
params = {'backend': 'ps',
          'font.size': 12,
          'axes.labelsize': 12,
         # 'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,    #benutze latex fuer schrift encoding -> automatisch selbe schriftart wie in latex
          'text.latex.unicode': True,
          'font.family': 'serif',
          'font.serif': 'cm',
          #'figure.figsize': fig_size,
          'text.latex.preamble': [r'\usepackage{physics}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)

plt.xlabel('$\mbox{dephasing} / s^{-1}$', fontsize = 60)
plt.ylabel('$P(\ket{-})$', fontsize = 60)
#plt.title('Measurement probabilities vs. dephasing', fontsize = 22)
#plt.grid(True)
plt.gcf().subplots_adjust(bottom=0.15, left = 0.15)
plt.legend()
plt.legend(loc='upper right', fontsize = 50)
plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1],[r'$0.4$',r'$0.5$', r'$0.6$', r'$0.7$', r'$0.8$', r'$0.9$', r'$1.0$', r'$1.1$', r'$1.2$'], fontsize = 60)
plt.tick_params(axis='both', which='major', pad=15)

plt.xticks([ -1, 0, 1, 2, 3, 4, 5, 6, 7, 8],[r'$10^{-1}$', r'$0$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$', r'$10^8$'], fontsize = 60)



ax = plt.gca()


""" Labels for each spin species
p0 = ax.text(circ_dephasing[0], 1, "$T^{e*}_2$, P (nat. Si, SET)", size=20, va="center", ha="center", rotation=30)
p1 = ax.text(circ_dephasing[1], 1, "$T^{e*}_2$, SiC (RT) ", size=20, va="center", ha="center", rotation=30)
p2 = ax.text(circ_dephasing[2], 1, "$T^{e*}_2, SiC (20 K)$", size=20, va="center", ha="center", rotation=30)
p3 = ax.text(circ_dephasing[3], 0.93, "$T^{e}_2, SiC (RT)$", size=20, va="center", ha="center", rotation=30)
p4 = ax.text(circ_dephasing[4], 0.95, "$T^{e*}_2$, P (puri. Si, SET)", size=20, va="center", ha="center", rotation=30)
p5 = ax.text(circ_dephasing[5], 1.0, "$T^{e}_2$, P (nat. Si, SET)", size=20, va="center", ha="center", rotation=30)
p6 = ax.text(circ_dephasing[6], 1.00, "$T^{e}_2$, P (puri. Si, SET)", size=20, va="center", ha="center", rotation=30)
p7 = ax.text(circ_dephasing[7], 1.05, "$T^{e}_2$, SiC (20 K)", size=20, va="center", ha="center", rotation=30)
p8 = ax.text(circ_dephasing[8], 1.1, "$T^{e}_2$, NV (puri. C, RT)", size=20, va="center", ha="center", rotation=30)
p9 = ax.text(circ_dephasing[9], 1.05, "$T^{e}_2$, Bi (puri. Si, CT)", size=20, va="center", ha="center", rotation=30)
"""

for i in range(0, 10):
    p0 = ax.text(circ_dephasing[i], abrupt_probs[i] + 0.04, '(' + str(10 -i) + ')', size=40, va="center", ha="center")

plt.grid(True, color='gray', linestyle='--', linewidth=2)
fig = plt.gcf()
st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')
plt.savefig('dephasing/phase_graph.png',  transparent=False, dpi=500)
plt.show()
