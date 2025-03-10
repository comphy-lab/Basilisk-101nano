/**
 * # 1D Transient Heat Conduction Solver (Basilisk version)
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
 */

#include "grid/cartesian1D.h"
#include "run.h"

// Declare scalar field for temperature
scalar T[];

// Simulation parameters
#define EPS 0.1  // Width of initial temperature peak
#define tmax 1.0
#define tsnap 0.1
double K;
int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  N = 200;       // Number of cells
  
  // Set the timestep based on stability criterion (CFL condition)
  // dt = dx^2/2 for explicit scheme
  K = XX; // fill in the value of K for a stable solution
  DT = (L0/N)*(L0/N)/K;
  
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
 * Time integration using explicit finite volume method
 */
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // TODO: Implement the heat equation integration
  // HINT: You need to compute fluxes at cell faces and update temperature field
  
  // 1. Define a scalar field for temperature fluxes
  scalar q[];
  
  // 2. Compute fluxes at cell faces
  // q_i = -(T_i - T_{i-1})/Delta
  // YOUR CODE HERE
  
  // 3. Compute temperature change rate
  // dT/dt = -(q_{i+1} - q_i)/Delta
  scalar dT[];
  // YOUR CODE HERE
  
  // 4. Update temperature field using explicit Euler step
  // T = T + dt * dT
  // YOUR CODE HERE
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
