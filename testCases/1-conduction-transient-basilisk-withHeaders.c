/**
 * # 1D Transient Heat Conduction Solver (Using Basilisk diffusion.h)
 * 
 * This program solves the transient heat conduction equation in one dimension:
 * 
 * \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}
 * 
 * Subject to no-flux boundary conditions on both ends of the domain.
 * Initial condition is a "Dirac delta" approximated by a thin rectangle
 * centered at x=0 with total integral = 1.
 * 
 * The exact self-similar analytical solution is:
 * T(x,t) = \frac{1}{2\sqrt{\pi t}}e^{-x^2/4t}
 * 
 * This version uses the diffusion.h module from Basilisk to handle the
 * diffusion equation implicitly.
 */

/* Include necessary headers in the correct order for Basilisk */
#include "grid/multigrid1D.h"  /* Multigrid solver is required by diffusion.h */
#include "run.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];
// Boundary conditions
// The diffusion solver will use homogeneous Neumann conditions by default
T[left] = neumann(0.);
T[right] = neumann(0.);

// Simulation parameters
#define EPS 0.1  // Width of initial temperature peak
#define tmax 1.0
#define tsnap 0.1

int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  N = 10000;       // Number of cells
  
  // We can use a larger timestep with the implicit solver
  // compared to the explicit version which requires dt = dx^2/2
  DT = 0.01;
  
  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

/**
 * Initialize temperature field
 * 
 * Sets up a "Dirac delta" approximated by a thin rectangle
 * centered at x=0 with total integral = 1.
 */
event init (t = 0) {
  foreach()
    T[] = (fabs(x) < EPS) ? 1.0/EPS/2.0 : 0.0; 
}

/**
 * Time integration using implicit diffusion solver
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Use the diffusion() function from diffusion.h to solve the equation
  // The heat equation is: ∂T/∂t = ∇²T which corresponds to diffusion with D = 1
  // face vector D[];
  // foreach_face()
  //   D.x[] = 1.0; // Constant diffusion coefficient of 1.0
  
  diffusion(T, dt);
}

/**
 * Save snapshots at regular intervals
 */
event writingFiles (t += tsnap; t < tmax+tsnap) {
  char filename[100];
  sprintf(filename, "intermediate/snapshot-%5.4f.csv", t);  

  FILE *fp = fopen(filename, "w");  
  foreach() {
    fprintf(fp, "%g,%g\n", x, T[]);
  }
  fclose(fp);
}

/**
 * Save final results and comparison with analytical solution
 */
event end (t = end) {
  char filename[100];
  sprintf(filename, "conduction-transient.csv");

  FILE *fp = fopen(filename, "w");
  foreach() {
    fprintf(fp, "%g,%g\n", x, T[]);
  }
  fclose(fp);
} 