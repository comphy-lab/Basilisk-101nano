/**
# Lid-Driven Cavity Flow of a Bingham Fluid

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
  /**
  ### Implementation of Generalized Newtonian Viscosity
  
  Rate-of-strain tensor components:
  - D₁₁ = ∂u/∂x
  - D₁₂ = D₂₁ = (1/2)(∂u/∂y + ∂v/∂x)
  - D₂₂ = ∂v/∂y
  
  Second invariant: D₂ = √(D_ij·D_ij)
  
  The equivalent viscosity is:
  $$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 }$$
  
  Note: ||D|| = D₂/√2
  
  We regularize the viscosity by limiting it to μ_max:
  $$\mu = \min(\mu_{eq}, \mu_{max})$$
  */
  foreach_face() {
    double D11 = (u.x[] - u.x[-1,0]); 
    double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
    double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + 
                      (u.y[] - u.y[-1,0]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    
    double mu_local;
    if (D2 > 0.) {
      double temp = tauy/(sqrt(2.0)*D2) + MU_0;
      mu_local = min(temp, mumax);
    } else {
      mu_local = (tauy > 0.0) ? mumax : MU_0;
    }
    
    muv.x[] = fm.x[]*(mu_local);
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

    fprintf(ferr, "Doing counter %d\n", counter);
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
      mumax = 1e2;
    }
    
    fprintf(ferr, "tauy = %g\n", tauy);
    sprintf(filename, "tau%d", counter);
    
    run();
  }
}