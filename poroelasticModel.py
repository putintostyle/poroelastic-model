import matplotlib.pyplot as plt
from numpy import select
from uniformGrid import *
from parameterSet import *
from simplify import *
uniformScheme = uniformGrid()
uniformScheme.IC()

simplifier = simplifier()
simplifier.parameters()
c1, c2, c3, c4 = simplifier.getParameter()
print(c1, c2, c3, c4)
u_steady = c3*uniformScheme.r_grid+c4/uniformScheme.r_grid**2+c1/2*simplifier.xi
p_steady = c1/uniformScheme.r_grid+c2
uniformScheme.u0 = u_steady
uniformScheme.p0 = p_steady

fig, (ax1, ax2, ax3) = plt.subplots(3)
# ax1.plot(uniformScheme.r_grid, uniformScheme.u0)
# ax2.plot(uniformScheme.r_grid, uniformScheme.u_dot)
# ax3.plot(uniformScheme.r_grid, uniformScheme.p0)
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
    
    ax1.plot(uniformScheme.r_grid, uniformScheme.u0-u_steady)
    ax2.cla()
    ax2.plot(uniformScheme.r_grid, uniformScheme.u_dot)
    ax3.cla()
    ax3.plot(uniformScheme.r_grid, uniformScheme.p0)
    plt.suptitle('t={}'.format(t*uniformScheme.dt))   
    
    plt.pause(uniformScheme.dt/100)
plt.show() 
