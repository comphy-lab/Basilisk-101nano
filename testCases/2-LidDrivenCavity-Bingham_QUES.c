/**
## Lid-Driven Cavity Flow of a Bingham Fluid

This simulation models a lid-driven cavity flow for a generalized Newtonian fluid,
specifically a Bingham fluid. The code uses regularization of the stress-strain
relationship instead of the Augmented Lagrangian Method.

The second invariant of the stress tensor is calculated at face centers.

## Parameters
- Newtonian: μ₀ = 1.0; τᵧ = 0.0 and n = 1
- Bingham: μ₀ = 1.0; τᵧ > 0.0 and n = 1
*/
#include "navier-stokes/centered.h"

// Constants
#define MU_0    (1.0)   // Base viscosity
#define LEVEL   6       // Grid refinement level
#define MAXDT   (1e-4)  // Maximum timestep
#define ERROR   (1e-6)  // Convergence error threshold

// Global variables
char filename[80];
double tauy;            // Yield stress
int counter;
double mumax = 1e3;     // Maximum regularized viscosity
int imax = 1e5;         // Maximum iterations

// Scalar and vector fields
scalar un[];            // Previous x-velocity (for convergence check)
face vector muv[];      // Face-centered viscosity field

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
  // Set viscosity field for non-Newtonian fluid
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
## Non-Newtonian Viscosity Properties
The viscosity is updated at each timestep based on the strain rate.
*/
event properties(i++) {
  foreach_face() {
    // TODO: Calculate deformation tensor components
    double D11 = /* Your code here */; 
    double D22 = /* Your code here */;
    double D12 = /* Your code here */;
    
    // TODO: Calculate second invariant (D2)
    double D2 = /* Your code here */;
    
    double mu_local;
    if (D2 > 0.) {
      // TODO: Implement regularized Bingham model
      double temp = /* Your code here */;
      mu_local = /* Your code here */;
    } else {
      // Handle zero strain rate case
      mu_local = (tauy > 0.0) ? mumax : MU_0;
    }
    
    // Apply viscosity
    muv.x[] = fm.x[] * (mu_local);
  }
}

/**
## Convergence Check
Monitor convergence every 500 iterations
*/
event logfile (i=i+500; i <= imax) {
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
  
  // Run simulations for different yield stress values
  for (counter = 0; counter < 5; counter++) {
    // Set yield stress for this run
    if (counter == 0) {
      tauy = 0.0;  // Newtonian case
    } else {
      tauy = pow(10, counter-1)/sqrt(2);
    }
    
    // Use higher viscosity cap for higher yield stress cases
    if (counter > 1) {
      mumax = 1e4;
    } else {
      mumax = 1e3;
    }
    
    fprintf(ferr, "tauy = %g\n", tauy);
    sprintf(filename, "tau%d", counter);
    
    // Store current velocity for convergence check
    foreach() {
      un[] = u.x[];
    }
    
    run();
  }
}