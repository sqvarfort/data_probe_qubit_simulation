#!/usr/bin/python
# -*- coding: latin-1 -*-
from qutip import *
import pickle
from pylab import *

import os
import sys
reload(sys)
sys.setdefaultencoding('utf-8')


#set image size for latex 
fig_width_pt = 276  # Get this from LaTeX using \showthe\linewidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]#in inch
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
          'figure.figsize': fig_size}
rcParams.update(params)


philist=pickle.load( open( "philist.pckl", "rb" ) )
taulist=pickle.load( open( "taulist.pckl", "rb" ) )

p=polyfit(taulist*1e3, philist, 1)

tau=(pi/2-p[1])/p[0]

fig=figure(1)
clf()
ax=fig.add_subplot(1,1,1)

ax.plot(taulist*1e3, philist, 'k.')
x=arange(3.e-3,3.6e-3,0.01e-3)*1e3
ax.plot(x, p[0]*x+p[1], 'r:')
ax.axhline(pi/2, color='0.1', ls=':')
ax.set_xlabel("complete circle time (ms)")
ax.set_ylabel("rotation (rad)")

subplots_adjust(left=0.15,right=0.95,top=0.90,bottom=0.2,hspace = .15,wspace=0.15)#hspace
savefig('circle_find_tau.svg')
show()
