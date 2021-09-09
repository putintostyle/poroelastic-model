import numpy as np
from parameterSet import *
class uniformGrid:
    def __init__(self):
        self.parameter = parameter
        self.M = 500
        self.T = 0.002
        self.dr = (self.parameter.r_s-self.parameter.r_v)/(self.M+1)
        self.dt = 1e-6
        self.r_grid = np.arange(self.parameter.r_v, self.parameter.r_s, self.dr )
        self.t_grid = np.arange(0, self.T, self.dt)

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
                   ((1-self.parameter.nu)((-3*self.u0[0]+4*self.u0[1]-self.u0[2])/(2*self.dr))+2*self.parameter.nu/self.parameter.r_v*self.u0[0])
        self.p0[0] = self.p_s
    
    def BC_dot(self):
        u_dot_s = 1/4*np.pi/(self.parameter.r_v+self.u0[0])**2*\
                        (self.parameter.Q_prod-np.pi*self.parameter.d**4/128/self.parameter.mu/self.parameter.L*\
                        (self.p0[0]-self.p0[-1])+4*np.pi*self.parameter.k*(self.parameter.r_v+self.u0[0])**2*(-3*self.p0[0]+4*self.p0[1]-self.p0[2])/(2*self.dr))
        
        self.u_dot[0] = u_dot_s
    def central_diff(self, array, diff_order, dx): # only in interior points, ie, if len(array) = N we only consider array[2::]-array[0::-2]
        if diff_order == 1:
            return (array[2::]-array[0:-2])/(2*dx)
        elif diff_order == 2:
            return (array[2::]-2*array[1:-1]+array[0:-2])/dx**2

    def update_u(self): # time
        self.u1 = self.u_dot*self.dt+self.u0 

    def update_u_dummy(self): # We only consider interior points
        self.u_dummy = np.zeros(len(self.u_dot))
        self.u_dummy[1:-1] = 2*self.parameter.G*(1-self.parameter.mu)/self.parameter.rho/(1-2*self.parameter.nu)\
                        *(self.central_diff(self.u1, 2, self.dr)
                            +self.central_diff(self.u1, 1, self.dr)/self.r_grid[1:-1]\
                            -(2/self.r_grid**2)*self.u1[1:-1])\
                        -self.parameter.alpha/self.parameter.rho*self.central_diff(self.p0, 1, self.dr)
        
        self.u_dummy[0] = None
        self.u_dummy[-1] = 0
    
    def update_p(self):
        p1 = np.zeros(len(self.p0))
        p1[1:-1] = self.p0[1:-1]+\
            (self.parameter.Q_prod)*\
            (self.parameter.alpha*(self.central_diff(self.u_dot, 1, self.dr)+2/self.dr*self.u_dot[1:-1])\
            +self.parameter.kappa*(self.central_diff(self.p0, 2, self.dr)+2/self.r_grid*self.central_diff(self.p0, 1, self.dr))\
            +self.parameter.kappa*self.parameter.rho(self.central_diff(self.u_dummy, 1, self.dr)+2/self.r_grid*self.u_dummy))
        
        p1[0] = 2*self.parameter.G/(1-2*self.parameter.nu)/(self.parameter.alpha-1)*\
                ((1-self.parameter.nu)*(-3*self.u0[0]+4*self.u0[1]-self.u0[2])/2/self.dr)+2*self.parameter.nu/self.r_grid[0]*self.u0[0]
        p1[-1] = self.p0[-1]
        self.p0 = p1