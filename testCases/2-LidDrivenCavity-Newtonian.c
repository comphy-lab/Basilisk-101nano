/**
# Lid-Driven Cavity Flow of a Newtonian Fluid

This simulation models a lid-driven cavity flow for a Newtonian fluid with
constant viscosity. This is a classic benchmark case in computational fluid dynamics.

## Parameters
- Reynolds number: Re = ρUL/μ = 1/μ (with ρ=1, U=1, L=1)
- We use μ = 1.0 by default (Re = 1)
*/
#include "navier-stokes/centered.h"

// Constants
#define LEVEL   8       // Grid refinement level
#define MAXDT   (1e-4)  // Maximum timestep
#define ERROR   (1e-6)  // Convergence error threshold

// Global variables
char filename[80] = "newtonian";  // Output filename
int imax = 1e5;                   // Maximum iterations

// Scalar field for convergence check
scalar un[];  // Previous x-velocity
const face vector muv[] = {1.0, 1.0};      // Face-centered viscosity field

/**
## Boundary Conditions
*/
// Top moving wall (lid)
u.t[top] = dirichlet(1);
// Other no-slip boundaries
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);
uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;

/**
## Initialization
*/
event init (t = 0) {
  // Set constant viscosity for Newtonian fluid
  mu = muv;
  
  // Initialize velocity field
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    un[] = 0;
  }
  
  dump (file = "start");
}

/**
## Convergence Check
Monitor convergence every 500 iterations
*/
event logfile (i=i+100; i <= imax) {
  double du = change(u.x, un);
  fprintf(ferr, "i = %d: dt = %g err = %g\n", i, dt, du);
  
  if (i > 0 && du < ERROR) {
    dump(file = filename);
    return 1; /* stop */
  }
  
  if (i > imax-10) {
    dump(file = filename);
  }
}

/**
## Output & Visualization
*/
event output (t = end) {  
  // Output fields in a format suitable for visualization
  dump(file="results");
}

/**
## Main Function
*/
int main() {
  // Initialize grid and parameters
  init_grid(1<<LEVEL);
  L0 = 1.0;
  origin(-0.5, -0.5);
  DT = MAXDT;
  TOLERANCE = 1e-5;
  CFL = 0.25;
  
  // Store current velocity for convergence check
  foreach() {
    un[] = u.x[];
  }
  
  // Run simulation
  run();
  
  return 0;
}
