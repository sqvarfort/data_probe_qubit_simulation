from qutip import *
from pylab import *
import pickle

philist=pickle.load( open( "philist2pi.pckl", "rb" ) )
taulist=pickle.load( open( "taulist2pi.pckl", "rb" ) )

p=polyfit(taulist, philist, 1)

#tau=7.7155254093495203e-05
#tau=0.00030859857755421354

tau=(2.*pi-p[1])/p[0]

print "2pi tau: %e" % tau

figure(1)
plot(taulist, philist, 'k.')
x=arange(4*72.e-6,4*83e-6,0.4e-6)
plot(x, p[0]*x+p[1], 'r:')
show()
