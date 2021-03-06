import numpy as np
from parameterSet import *

class uniformGrid:
    def __init__(self):
        self.parameter = parameter()
        self.M = self.parameter.M
        self.T = self.parameter.T
        self.dr = (self.parameter.r_s-self.parameter.r_v)/(self.M+1)
        self.dt = self.parameter.dt
        self.r_grid = np.linspace(self.parameter.r_v, self.parameter.r_s, self.M+1)
        self.t_grid = np.arange(0, self.T, self.dt)
        self.u0 = None
        self.p0 = None
        self.u_dummy_coefficient_1 = 2*self.parameter.G*(1-self.parameter.nu)/self.parameter.rho/(1-2*self.parameter.nu)
        self.u_dummy_coefficient_2 = self.parameter.alpha/self.parameter.rho

        self.u_dot_coefficient = self.parameter.d**4*np.pi/128/self.parameter.mu/self.parameter.L
    def IC(self, impact = True):
        
        if impact:
            # self.impactFunction = impactFunction()
            k = np.floor(self.M/3.5)
            amp = 0.004
            self.u_dot = amp*np.sin(100*self.r_grid*np.pi)
            
            self.u_dot[:int((self.M-np.floor(k/2)))] = 0
        else:
            self.u_dot = np.zeros(self.M+1)
        # print(len(self.u_dot))


    def central_diff(self, array, diff_order, dx): # only in interior points, ie, if len(array) = N we only consider array[2::]-array[0::-2]
        if diff_order == 1:
            return (array[2::]-array[0:-2])/(2*dx)
        elif diff_order == 2:
            return (array[2::]-2*array[1:-1]+array[0:-2])/(dx**2)

    def update_u(self): # time
        self.u1 = self.u_dot*self.dt+self.u0 
        self.u1[-1] = 0
        

    def update_u_dummy(self): # We only consider interior points why consider u1
        self.u_dummy = np.zeros(self.M+1)
        self.u_dummy[1:-1] = self.u_dummy_coefficient_1\
                            *(self.central_diff(self.u0*self.r_grid, 2, self.dr)/self.r_grid[1:-1]\
                            -(2/(self.r_grid[1:-1]**2))*self.u0[1:-1])\
                            -self.u_dummy_coefficient_2*self.central_diff(self.p0, 1, self.dr)
        
        self.u_dummy[0] = self.u_dummy_coefficient_1\
                          *(np.sum(np.array([-2, 5, -4, 1])*self.r_grid[0:4]*self.u0[0:4])/(-1*self.dr**2)/self.r_grid[0]\
                          -2/(self.r_grid[0]**2)*self.u0[0])\
                          -self.u_dummy_coefficient_2*(-3*self.p0[0]+4*self.p0[1]-self.p0[2])/2/self.dr
        
        self.u_dummy[-1] = 0
        
        
    def update_p(self):
        
        p1 = np.zeros(self.M+1)
        p1[1:-1] = self.p0[1:-1]+\
            self.dt/(self.parameter.S)*\
            (-1*self.parameter.alpha*(self.central_diff(self.u_dot*self.r_grid**2, 1, self.dr)/self.r_grid[1:-1])\
            +self.parameter.k*(self.central_diff(self.p0*self.r_grid, 2, self.dr)/self.r_grid[1:-1])\
            +self.parameter.k*self.parameter.rho*(self.central_diff(self.u_dummy*(self.r_grid)**2, 1, self.dr)/self.r_grid[1:-1]**2))
        
        p1[0] = 2*self.parameter.G/(1-2*self.parameter.nu)/(self.parameter.alpha-1)*\
                ((1-self.parameter.nu)*(-3*self.u1[0]+4*self.u1[1]-self.u1[2])/2/self.dr+2*self.parameter.nu/self.r_grid[0]*self.u1[0])
        
        p1[-1] = self.parameter.p_v+self.parameter.mu*self.parameter.R*self.parameter.Q_obs
        self.p1 = p1
        
    def update_u_dot(self):
        
        u_dot_1 = np.zeros(self.M+1)
        u_dot_1[1:-1] = self.u_dot[1:-1]+\
                        self.dt*(self.u_dummy_coefficient_1\
                        *(self.central_diff(self.u1*self.r_grid, 2, self.dr)/self.r_grid[1:-1]-2/(self.r_grid[1:-1]**2)*self.u1[1:-1])\
                        -self.u_dummy_coefficient_2*self.central_diff(self.p1, 1, self.dr)
                        )
        u_dot_1[0] = 1/4/(np.pi*(self.r_grid[0]+self.u1[0])**2)\
                        *((self.parameter.Q_prod-(self.u_dot_coefficient)*(self.p1[0]-self.p1[-1])\
                            +4*self.parameter.k*np.pi*(self.r_grid[0]+self.u1[0])**2\
                            *(-3*self.p1[0]+4*self.p1[1]-self.p1[2])/2/self.dr))
        
        u_dot_1[-1] = 0
        self.u_dot_1 = u_dot_1

    