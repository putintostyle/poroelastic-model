# Current Problem
## Equations:
$$
\left\{
    \begin{array}{lll}
        \displaystyle\frac{\partial u}{\partial t} & = & \dot{u} \\[10pt]
        \displaystyle \frac{\partial \dot{u}}{\partial t} & = & \displaystyle \frac{2G(1-\nu)}{\rho(1-2\nu)}\left(\frac{\partial ^2 u}{\partial  r^2}+\frac{2}{r}\frac{\partial  u}{\partial  r}-\frac{2}{r^2}u\right)-\frac{\alpha}{\rho}\frac{\partial p}{\partial r}\\[10pt]
        \displaystyle
        \frac{\partial p}{\partial t} & = & \displaystyle Q^e\left(-
        \alpha^e \left(\frac{\partial  \dot{u}}{\partial  r}+\frac{2}{r}\dot{u}\right)+\kappa^e\left(\frac{\partial^2 p}{\partial^2 r} +\frac{2}{r}\frac{\partial p}{\partial r}\right)+\kappa^e\rho^e\left(\frac{\partial }{\partial r}\left(\frac{\partial \dot{u}}{\partial t}\right)+\frac{2}{r}\frac{\partial \dot{u}}{\partial t}\right)\right)
    \end{array}
\right.
$$

## Numerical Scheme:

Consider the first time step, wthat is $t = \Delta t$
* Update $u$ first:
    $$
        \displaystyle \frac{u^{1}_j-u^0_j}{\Delta t} = \dot{u}^0_j.   
    $$
* Update $p$ next:
    $$
        \displaystyle\frac{1}{Q^e}\frac{p^{1}_j-p^0_j}{\Delta t} =  \displaystyle-
        \alpha^e\left( \frac{\dot{u}^0_{j+1}-\dot{u}^0_{j-1}}{2\Delta r}+\frac{2}{r_j}\dot{u}^0_j\right)+\kappa^e\left(\frac{p^0_{j+1}-2p^0_{j}+p^0_{j-1}}{(\Delta r)^2} +\frac{2}{r_j}\frac{p^0_{j+1}-p^0_{j-1}}{\Delta r}\right)+\kappa\rho\frac{\left(\frac{\partial\dot{u}}{\partial t}\right)^0_{j+1}-\left(\frac{\partial\dot{u}}{\partial t}\right)^0_{j-1}}{2\Delta r}+\frac{2}{r_j}\left(\frac{\partial\dot{u}}{\partial t}\right)^0_{j}
    $$
    In addition, we have
    $$
        \left(\frac{\partial\dot{u}}{\partial t}\right)^0_{j} = \displaystyle\frac{2G(1-\nu)}{\rho(1-2\nu)}\left(\frac{u^{0}_{j+1}-2u^{0}_{j}+u^{0}_{j-1}}{(\Delta r)^2} +\frac{2}{r_j}\frac{u^{0}_{j+1}-u^{0}_{j-1}}{2\Delta r}-\frac{2}{r_j^2}u^{0}_j\right)-\frac{\alpha}{\rho}\frac{p^{0}_{j+1}-p^{0}_{j-1}}{2\Delta r}
    $$
    ```Question:```

    * 要更新 $p_1$，我們需要$\left(\frac{\partial\dot{u}}{\partial t}\right)^0_{0}$，但是我們求不出來（管哥直接用外插去算在邊界的微分）


* Finally updata $\dot{u}$:
    $$
        \displaystyle \frac{\dot{u}^{1}_j-\dot{u}^0_j}{\Delta t}  =  \displaystyle\frac{2G(1-\nu)}{\rho(1-2\nu)}\left(\frac{u^{1}_{j+1}-2u^{1}_{j}+u^{1}_{j-1}}{(\Delta r)^2} +\frac{2}{r_j}\frac{u^{1}_{j+1}-u^{1}_{j-1}}{2\Delta r}-\frac{2}{r_j^2}u^{1}_j\right)-\frac{\alpha}{\rho}\frac{p^{1}_{j+1}-p^{1}_{j-1}}{2\Delta r}
    $$