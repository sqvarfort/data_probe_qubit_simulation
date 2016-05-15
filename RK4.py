import numpy as np
from numpy import exp


def RK4(Lindblad, rho, h, time):
    k1 = Lindblad(time, rho)
    k2 = Lindblad(time + h/2., rho + h/2.*k1)
    k3 = Lindblad(time + h/2., rho + h/2.*k2)
    k4 = Lindblad(time + h, rho + h*k3)
    return (k1 + 2.*k2 + 2.*k3 + k4)*h/6.0



sigmax = np.matrix([[0,1], [1,0]]) # This is our

rho = np.matrix([[1,1],[1,1]])*1/2.


def commutator(A,B):
    return np.dot(A,B) - np.dot(B,A)

def Heisenberg(time, rho):
    sigmaz = np.matrix([[1,0],[0,-1]])
    return -1j*commutator(sigmaz,rho)



#def rho(t):
#    return np.matrix([[1,exp(-1j*t)],[exp(1j*t),1]])*1/2

h = 0.1
time = 0.0


import matplotlib.pyplot as plt
plt.plot()
plt.ylabel('some numbers')


for t in range(0,20):
    dy = RK4(Heisenberg, rho, h, time)
    rho = rho + dy
    time = time + h
    print rho
    plt.scatter(time, rho[0,1])
    plt.scatter(time, rho[1,0], color = 'red')


plt.show()
