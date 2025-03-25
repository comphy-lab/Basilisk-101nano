/**
 # 2D Transient Heat Conduction Solver
 
 This program solves the transient heat conduction equation in two dimensions:
 
 $$
 \frac{\partial T}{\partial t} = \nabla^2 T
 $$
 
 ## Subject to boundary conditions:

  - T = 1 on the top boundary
  - T = 0 on the bottom, left, and right boundaries
 
 ## Initial condition:
  - T = 0 everywhere
 
 ## Question:
 
 Fill in this code. For answer, see [1-conduction-2D.c](1-conduction-2D.c)

*/

#include "run.h"
#include "diffusion.h"