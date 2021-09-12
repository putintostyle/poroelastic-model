import matplotlib.pyplot as plt
from numpy import select
from uniformGrid import *
from parameterSet import *

setting = parameter()
uniformScheme = uniformGrid()
## set parameters
setting.M = 500
setting.T = 0.002
setting.dt = 1e-6

uniformScheme.IC()


for t in range(int(setting.T/setting.dt)):
    uniformScheme.update_u()
    uniformScheme.update_u_dummy()
    uniformScheme.update_p()
    uniformScheme.update_u_dot()

    uniformScheme.u0 = uniformScheme.u1
    uniformScheme.u_dot = uniformScheme.u_dot_1
    uniformScheme.p0 = uniformScheme.p1
    

