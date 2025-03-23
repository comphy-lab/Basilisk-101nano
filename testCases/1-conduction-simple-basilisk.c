/**
 ## 1D Steady Heat Conduction Solver (Basilisk version)
 
 This program solves the heat conduction equation to reach steady state:
 
 $$
 \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}
 $$
 
 ## Subject to Dirichlet boundary conditions:
 
 - T(0) = 0
 - T(1) = 1
 
 The steady-state analytical solution is the linear profile $T(x) = x$.
*/

#include "grid/cartesian1D.h" // 1D grid
#include "run.h" // utility functions for running the simulation

// Declare scalar field for temperature
scalar T[];
T[left] = dirichlet(0.0); // left boundary
T[right] = dirichlet(1.0); // right boundary

int main() {
  // Domain setup
  L0 = 1.0;     // Domain length
  X0 = 0.0;    // Left boundary
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
 
 Starting with zero temperature everywhere, letting the boundary
 conditions drive the evolution toward steady state.
*/
event init (t = 0) {
  foreach()
    T[] = 0.0; 
}

/**
 ## Time integration using explicit finite volume method
 */
event marching (i++) {
  foreach() {
    // Proper heat equation time stepping
    T[] += DT*(T[1] - 2*T[] + T[-1])/(Delta*Delta);
  }
}

/**
 ## Monitor convergence to steady state
 */
event testingConvergence (i++; i < 100000) {
  double error = 0.0;
  foreach() {
    // Calculate error relative to analytical solution T(x) = x
    error += fabs(T[] - x);
  }
  error = error/N; // Normalize by number of cells
  
  // Only print error every 1000 iterations to avoid flooding console
  if (i % 1000 == 0) {
    printf("i: %d, Error: %g\n", i, error);
  }
  
  // Stop if we've reached steady state
  if (error < 1e-10) {
    printf("Converged at iteration %d with error %g\n", i, error);
    return 1; // Trigger end of simulation
  }
}

/**
 ## Save final results and comparison with analytical solution
 */
event end (t = end) {
  char filename[100];
  sprintf(filename, "conduction-simple.csv");

  FILE *fp = fopen(filename, "w");
  foreach() {
    fprintf(fp, "%g,%g\n", x, T[]);
  }
  fclose(fp);
}
