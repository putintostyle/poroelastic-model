import matplotlib.pyplot as plt
from numpy import select
from uniformGrid import *
from parameterSet import *

uniformScheme = uniformGrid()
uniformScheme.IC()

fig, (ax1, ax2, ax3) = plt.subplots(3)

ax1.plot(uniformScheme.r_grid, uniformScheme.u0)
ax2.plot(uniformScheme.r_grid, uniformScheme.u_dot)
ax3.plot(uniformScheme.r_grid, uniformScheme.p0)
# plt.show()

for t in range(600):
    
    uniformScheme.update_u()
    uniformScheme.update_u_dummy()
    uniformScheme.update_p()
    uniformScheme.update_u_dot()

    uniformScheme.u0 = uniformScheme.u1
    uniformScheme.u_dot = uniformScheme.u_dot_1
    uniformScheme.p0 = uniformScheme.p1
    
    ax1.cla()
    
    ax1.plot(uniformScheme.r_grid, uniformScheme.u0)
    ax2.cla()
    ax2.plot(uniformScheme.r_grid, uniformScheme.u_dot)
    ax3.cla()
    ax3.plot(uniformScheme.r_grid, uniformScheme.p0)
    plt.suptitle('t={}'.format(t*uniformScheme.dt))   
    
    plt.pause(0.01)
plt.show() 
