/**
 ## 1D Transient Heat Conduction Solver (Basilisk version)
 
 This program solves the transient heat conduction equation in one dimension:
 
 $$
 \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}
 $$
 
 - Subject to no-flux boundary conditions on both ends of the domain.
 - Initial condition is a "Dirac delta" approximated by a thin rectangle centered at $x=0$ with total integral = 1.
 
 The exact self-similar analytical solution is:
 
 $$
 T(x,t) = \frac{1}{2\sqrt{\pi t}}e^{-\frac{x^2}{4t}}
 $$
*/

#include "grid/cartesian1D.h"
#include "run.h"

// Declare scalar field for temperature
scalar T[];

// Simulation parameters
#define EPS 0.1  // Width of initial temperature peak
#define tmax 1.0
#define tsnap 0.1

int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  N = 200;       // Number of cells
  
  // Set the timestep based on stability criterion (CFL condition)
  // dt = dx^2/2 for explicit scheme
  DT = (L0/N)*(L0/N)/2;
  
  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

/**
 ## Initialize temperature field
 
 Sets up a `Dirac delta` approximated by a thin rectangle centered at $x=0$ with total integral = 1.
*/
event init (t = 0) {
  foreach()
    T[] = (fabs(x) < EPS) ? 1.0/EPS/2.0 : 0.0; 
}

/**
 ## Time integration using explicit finite volume method
*/
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Compute fluxes at the faces of the cells
  // q_i = -(T_i - T_{i-1})/Delta
  scalar q[];
  foreach()
    q[] = -(T[] - T[-1])/Delta;
  
  // Update temperature
  // dT/dt = -(q_{i+1} - q_i)/Delta
  // this enforces the central difference scheme
  scalar dT[];
  foreach()
    dT[] = -(q[1] - q[])/Delta;
  
  // Explicit Euler step
  foreach()
    T[] += dt*dT[];
  
}

/**
 ## Save snapshots at regular intervals
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
 ## Save final results and comparison with analytical solution
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
