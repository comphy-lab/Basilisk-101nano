/**
#  Planar Couette flow of Newtonian Fluid
# Code
*/
#include "navier-stokes/centered.h"

// Global parameters
char file_name[80];
double tau_y = 0.0;        // Yield stress
double mu_0 = 1.0;         // Base viscosity
int max_iter = 1e4;        // Maximum iterations
#define DT_MAX (1e-3)      // Maximum timestep

int main(int argc, char const *argv[])
{
  
  sprintf(file_name, "results");
  
  // Initialize grid and domain
  init_grid(1<<8);
  L0 = 1.0;
  origin(-0.5, -0.5);
  DT = DT_MAX;
  stokes = true;
  TOLERANCE = 1e-5;

/**
 Values of yield stress, viscosity, and coefficient.
 - Newtonian: $\mu_0 = 1.0$
*/
  // Set boundary conditions
  // Right-left boundaries are periodic
  periodic(right);
  
  // Slip at the top
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  uf.n[top] = neumann(0);
  
  // No slip at the bottom
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  
  // Pressure conditions are neumann 0.0
  p[top] = neumann(0);
  pf[top] = neumann(0);
  p[bottom] = neumann(0);
  pf[bottom] = neumann(0);

  run();
}

// un is used to search for a stationary solution
scalar un[];
// muv will be used as the face vector for viscosity
face vector muv[];

/**
## Initialization event
*/
event init(t = 0) {
  // Set Non-Newtonian viscosity
  mu = muv;
  
  /**
   Pressure gradient `mdpdx`
   $$-\frac{dp}{dx} = 1 $$
  */
  const face vector mdpdx[] = {1.0, 0.0};
  a = mdpdx;
  
  // Initialize velocity field at rest
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    un[] = 0;
  }
  
  dump(file = "start");
}

/**
## Monitoring convergence 
We look for a stationary solution by checking changes in velocity field.
*/
event logfile(i += 500; i <= max_iter) {
  double du = change(u.x, un);
  fprintf(ferr, "i = %d: err = %g\n", i, du);
  
  if (i > 0 && du < 1e-6) {
    dump(file = file_name);
    return 1; /* stop */
  }
  
  if (i == max_iter) {
    dump(file = file_name);
  }
}

/**
## Calculating viscosity for generalized Newtonian fluid
*/
event properties(i++) {
  foreach_face() {
    // Apply viscosity at face
    muv.x[] = fm.x[] * mu_0;
  }
}