from qutip import *
from pylab import *
import pickle

philist=pickle.load( open( "philist.pckl", "rb" ) )
taulist=pickle.load( open( "taulist.pckl", "rb" ) )

p=polyfit(taulist, philist, 1)

tau=(pi/2.-p[1])/p[0]

figure(1)
plot(taulist, philist, 'k.')
x=arange(3.e-3,3.6e-3,0.01e-3)
plot(x, p[0]*x+p[1], 'r:')
show()