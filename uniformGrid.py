import numpy as np
from parameterSet import *
class uniformGrid:
    def __init__(self):
        self.parameter = parameter
        self.M = 500
        self.T = 0.002
        self.dx = (self.parameter.r_s-self.parameter.r_v)/(self.M+1)
        self.dt = 1e-6
        
    def IC(self):
        self.u0 = np.zeros(self.M+1)
        self.p0 = (self.parameter.p_v+self.parameter.mu*self.parameter.R*self.parameter.Q_obs)*np.ones(self.M+1)
        self.u_dot = np.zeros(self.M+1)
    
    def BC(self):
        self.u_e = 0

        self.u_dot_e = 0

        self.p_e = self.p0

    def BC_p(self):
        self.p_s = (2*self.parameter.G)/(1-2*self.parameter.nu)/(self.parameter.alpha-1)\
                   ((1-self.parameter.nu)((-3*self.u0[0]+4*self.u0[1]-self.u0[2])/(2*self.dx))+2*self.parameter.nu/self.parameter.r_v*self.u0[0])
        self.p0[0] = self.p_s
    
    def BC_dot(self):
        self.u_dot_s = 1/4*np.pi/(self.parameter.r_v+self.u0[0])**2*\
                        (self.parameter.Q_prod-np.pi*self.parameter.d**4/128/self.parameter.mu/self.parameter.L*\
                        (self.p0[0]-self.p0[-1])+4*np.pi*self.parameter.k*(self.parameter.r_v+self.u0[0])**2*(-3*self.p0[0]+4*self.p0[1]-self.p0[2])/(2*self.dx))
        self.dot[0] = self.u_dot_s
    
    def update_u(self):
        self.u0 = self.u_dot*self.dt+self.u0
    
    def update_u_dummy(self):
        self.u_dummy = 

