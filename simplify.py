from parameterSet import *
from sympy import *
import numpy as np
class simplifier:
    def __init__(self):
        # self.Q_prod, self.d, self.mu, self.L, self.k = symbols('Q_prod, d, mu, L, kappa')
        self.c1, self.c2, self.c3, self.c4 = symbols('c1 c2 c3 c4')

        self.parameter = parameter()
        self.r_v = self.parameter.r_v
        self.r_s = self.parameter.r_s
        self.nu = self.parameter.nu
        self.alpha = self.parameter.alpha
        self.xi = ((1-2*self.parameter.nu)*self.parameter.alpha)/(2*self.parameter.G*(1-self.parameter.nu))

        self.p_s = self.parameter.p_v+self.parameter.mu*self.parameter.R*self.parameter.Q_obs
    def parameters(self):
        c1, c2, c3, c4, xi, r_v, r_s, nu, alpha, p_s = self.c1, self.c2, self.c3, self.c4, self.xi, self.r_v, self.r_s, self.nu, self.alpha, self.p_s
        # c1 = r_s*p_s-r_s*c2
        c2 = (r_s*p_s-c1)/r_s
        c4 = -(c3*r_s**3+(c1*xi*r_s**2)/2)
        eq4 = c3-2*c4/(r_v**3)+nu/(1-nu)*2/r_v*(c3*r_v+c4/r_v**2+c1*xi/2)-(1-1/alpha)*(c1/r_v+c2)
        c3 = solve(eq4, c3)[0]

        print(simplify(-(c3*r_s**3+(c1*xi*r_s**2)/2)))
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = -(c3*r_s**3+(c1*xi*r_s**2)/2)

    def solve_c1(self):
        
        Q_prod, d, mu, L, k = self.parameter.Q_prod, self.parameter.d, self.parameter.mu, self.parameter.L, self.parameter.k
        equation = Q_prod-(np.pi*d**4/128/mu/L*(self.c1/self.r_v+self.c2-self.p_s)\
                   -4*np.pi*k*(self.r_v+self.c3*self.r_v+self.c4/self.r_v**2+self.c1*self.xi/2)**2*(-1*self.c1/self.r_v**2))
        self.c1 = solve(equation, self.c1)[0]
        

    def getParameter(self):
        self.solve_c1()
        self.c2 = (self.parameter.r_s*self.p_s-self.c1)/self.parameter.r_s
        c3, c4 = symbols('c3 c4')
        xi, r_v, r_s, nu, alpha, p_s = self.xi, self.r_v, self.r_s, self.nu, self.alpha, self.p_s
        c1 = self.c1
        # c1 = r_s*p_s-r_s*c2
        c2 = (r_s*p_s-c1)/r_s
        c4 = -(c3*r_s**3+(c1*xi*r_s**2)/2)
        eq4 = c3-2*c4/(r_v**3)+nu/(1-nu)*2/r_v*(c3*r_v+c4/r_v**2+c1*xi/2)-(1-1/alpha)*(c1/r_v+c2)
        
        self.c3 = solve(eq4, c3)[0]
        self.c4 = -(self.c3*r_s**3+(c1*xi*r_s**2)/2)
        return self.c1, self.c2, self.c3, self.c4
# print(simplify(-alpha*c1*nu*r_s*xi - alpha*c1*nu*r_s + alpha*c1*nu*r_v*xi/2 + alpha*c1*nu*r_v + alpha*c1*r_s + alpha*c1*r_v*xi/2 - alpha*c1*r_v - alpha*nu*p_s*r_s*r_v + alpha*p_s*r_s*r_v + c1*nu*r_s - c1*nu*r_v - c1*r_s + c1*r_v + nu*p_s*r_s*r_v - p_s*r_s*r_v))
