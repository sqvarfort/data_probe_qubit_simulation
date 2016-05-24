from scipy import constants as cp
import numpy as np
from qutip import *
import yaml

class H_RWA(object):
    def __init__(self, config):
        """
            config can be file -> yaml or dict
            config need to be of form
            {'J': J, 'Delta': Delta, 'r': array([0,0,d]), 'circ': False, 'cOpts': {'pJit': False, 'rstd': 1.e-9, 'cOffset': zeros(3), 'tau': tau, 'd': d, 'D': D} }
            TODO: check with yaml list vs numpy array!
        """
        # needed time independent diagonal and off diagonal parts. diagonal part is just sigmazz. off diagonal part is H_12 and H_21 part of the full matrix.
        self.args={'Hd': tensor(sigmaz(), sigmaz()), 'H12': Qobj([[0,0,0,0],[0,0,0,0],[0,1,0,0],[0,0,0,0]]), 'H21': Qobj([[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]])}

        # load config
        if type(config) == str:
            with open(config) as cfile:
                self.args.update(yaml.load(cfile))
        elif type(config) == dict:
            self.args.update(config)
        else:
            print "config is wrong type"

        # Zeeman aprt in rotating frame of probe qubit. Probe qubit is the first qubit so only the sigmaz on data qubit remains in the zeeman hamiltonian
        sigma2z = tensor(qeye(2), sigmaz())
        self.args['Hz']= self.args['Delta']*sigma2z

    @staticmethod
    def circ_motion(t, args):
        # cOffset can be generated using Gavins function!
        r=np.array(args['cOffset']) + np.array( [ -args['D']/2. + args['D']/sqrt(2)*np.sin(2*np.pi/args['tau'] * t)  ,  -args['D']/2. + args['D']/np.sqrt(2)*cos(2*np.pi/args['tau'] * t)  ,  args['d'] ] )
        #allow random path jitter (simulates not perfect movement of mems stage) ~~1nm?
        if args['pJit']:
            r+=np.random.normal(0, args['rstd'], 3)
            # use below for test plot
            #r=args['cOffset'] + array( [ -args['D']/2. + random.normal(args['D']/sqrt(2)*sin(2*pi/args['tau'] * t), args['rstd'])  ,  -args['D']/2. + random.normal(args['D']/sqrt(2)*cos(2*pi/args['tau'] * t), args['rstd'])  ,  random.normal(args['d'], args['rstd']) ] )
        return r

    @staticmethod
    def getDelta(mat1='Bi', mat2='P', Bfield=300e-3):
        """
            mat1 anod mat2 is string of dict index
            return Delta in units of frequency
        """

        muB=cp.e*cp.hbar/(2*cp.m_e) #9.27e-24
        g={}
        g['e']=1.9985
        g['P']=g['e']-2.5e-4
        g['Bi']=2.0003

        return (g[mat2]-g[mat1]) * muB * Bfield / cp.hbar

    @staticmethod
    def getDeltaA(mat1, mat2):
        """
            mat1 mat 2 is string of dict index
            returns Delta inunits of freq calculated using hyperfine coupling strength
        """
        # account for nuclear spin -> 9/4, 1/4
        A={'P': 118.e6 / 4, 'Bi': 9*1475.4e6 / 4}

        return A[mat1] - A[mat2]

    @staticmethod
    def getJ():
        """
            returns J in units of freq
        """
        ge=1.9985
        muB=cp.e*cp.hbar/(2*cp.m_e) #9.27e-24
        J=cp.mu_0*ge**2 * muB**2/(4*np.pi)/cp.hbar #in units of frequency
        return J

    @staticmethod
    def __Hfunc(t, args):
        """
        ## Arguments need to be of form:
        # needed time independent diagonal and off diagonal parts. diagonal part is just sigmazz. off diagonal part is H_12 and H_21 part of the full matrix.
        args={'Hz': Hz, 'Hd': tensor(sigmaz(), sigmaz()), 'H12': Qobj([[0,0,0,0],[0,0,0,0],[0,1,0,0],[0,0,0,0]]), 'H21': Qobj([[0,0,0,0],[0,0,1,0],[0,0,0,0],[0,0,0,0]])}

        # specify simulation arguments
        args.update({'J': J, 'Delta': Delta, 'r': array([0,0,d]), 'circ': False, 'cOpts': {'pJit': False, 'rstd': 1.e-9, 'cOffset': zeros(3), 'tau': tau, 'd': d, 'D': D} })
        """
        r=np.array(args['r'])
        phase=np.exp(4j*args['Delta'])
        if args['circ']:
            #implements circular motion from with radius D/sqrt2 where D is the separation between donors. data qubit is in the origin with total circle time of tau!
            # cOffset allows to move not perfectly above the data qubit
            r=circ_motion(t, args['cOpts'])
        return args['Hz'] + args['J']/np.linalg.norm(r)**3 * ( 1. -3. * r[2]**2/np.linalg.norm(r)**2 ) * args['Hd'] + args['J']/np.linalg.norm(r)**3 * ( 2. -3.*r[0]**2/np.linalg.norm(r)**2 -3.*r[1]**2/np.linalg.norm(r)**2 ) * (phase*args['H12'] + np.conjugate(phase)*args['H21'])

    def getHfunc(self):
        return self.__Hfunc

    def getArgs(self):
        return self.args
