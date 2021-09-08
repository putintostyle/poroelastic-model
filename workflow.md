# ToDo List
### Implement codes
* Uniform Grid:
  * Task 1
    * space: Central difference
    * time : Backward Euler
    * 
  * Task 2
    * space : Central
    * time : Crank-Nicholsan
  * Check stability
* Dual Grid:
  * Do the same tasks as uniform grid has done
  
### Algroithm:
* Initial $u^0, p^0, v^0 = \dot{u}^0$
* Update $u^1$, update $p^1$ from $u^0$ $v^0$
  * We have $p^1(r_1) = \displaystyle\frac{2G}{(1-2\nu)(\alpha-1)}\left((1-nu)\left(\frac{-3u^1_0+4u^1_1-u^1_2}{2\Delta r}\right)+\frac{2\nu}{r_0}u^1_0\right)$
  * interior points: we first compute $v^*$ from $u^1$ and $p^0$
      * boundary : 外插
* Update $v^1$ from $u^1, p^1$