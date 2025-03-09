#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

/**
 * # 1D Transient Heat Conduction Solver
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

// Simulation parameters
#define N       200     // Number of cells
#define L0      10.0    // Domain length (from -L0/2 to +L0/2)
#define X0      (-L0/2) // Left boundary

/**
 * Create directory for intermediate output files
 * 
 * @return 0 on success, 1 on failure
 */
int create_output_directory() {
  struct stat st = {0};
  
  if (stat("intermediate", &st) == -1) {
    if (mkdir("intermediate", 0700) == -1) {
      fprintf(stderr, "Error creating intermediate directory: %s\n", strerror(errno));
      return 1;
    }
  }
  
  return 0;
}

/**
 * Initialize temperature array with a "Dirac delta" approximation
 * using a thin rectangle centered at x=0 with total integral = 1
 * 
 * @param temperature The temperature array to initialize
 * @param dx Cell size
 * @param eps Half-width of the initial rectangle
 */
void initialize_temperature(double *temperature, double dx, double eps) {
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx; // Cell-center coordinate
    if (fabs(x) < eps) {
      temperature[i] = 1.0 / (2.0 * eps);
    } else {
      temperature[i] = 0.0;
    }
  }
}

/**
 * Print current temperature field to console
 * 
 * @param temperature The current temperature array
 * @param dx Cell size
 * @param time Current simulation time
 */
void print_temperature(double *temperature, double dx, double time) {
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx;
    printf("%g  %g  %g\n", x, temperature[i], time);
  }
  printf("\n\n");
}

/**
 * Save a snapshot of the temperature field to a CSV file
 * 
 * @param temperature The current temperature array
 * @param dx Cell size
 * @param time Current simulation time
 * @return 0 on success, 1 on failure
 */
int save_snapshot(double *temperature, double dx, double time) {
  char filename[100];
  sprintf(filename, "intermediate/snapshot-%5.4f.csv", time);
  
  FILE *snap_file = fopen(filename, "w");
  if (snap_file == NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return 1;
  }
  
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx;
    fprintf(snap_file, "%g,%g\n", x, temperature[i]);
  }
  
  fclose(snap_file);
  printf("Saved snapshot at time = %5.4f\n", time);
  
  return 0;
}

/**
 * Compute heat fluxes at cell interfaces
 * q[i+0.5] = - (T[i+1] - T[i])/dx with no-flux boundary conditions
 * 
 * @param temperature Current temperature array
 * @param flux Array to store computed fluxes (size N+1)
 * @param dx Cell size
 */
void compute_fluxes(double *temperature, double *flux, double dx) {
  flux[0] = 0.0; // Left boundary (no flux)
  flux[N] = 0.0; // Right boundary (no flux)
  
  for (int i = 0; i < N - 1; i++) {
    flux[i+1] = - (temperature[i+1] - temperature[i]) / dx;
  }
}

/**
 * Update temperature field based on fluxes
 * 
 * @param t_current Current temperature array
 * @param t_new New temperature array to be computed
 * @param flux The computed flux array
 * @param dx Cell size
 * @param dt Current time step
 */
void update_temperature(double *t_current, double *t_new, double *flux, 
                        double dx, double dt) {
  for (int i = 0; i < N; i++) {
    double dq;
    
    if (i == 0) {
      // Leftmost cell: flux difference is q[0.5] - q[-0.5], but q[-0.5]=0
      dq = flux[1];
    } else if (i == N-1) {
      // Rightmost cell: flux difference is q[N-0.5] - q[N-1.5], but q[N+0.5]=0
      dq = - flux[N-1];
    } else {
      dq = flux[i+1] - flux[i];
    }
    
    t_new[i] = t_current[i] - (dt/dx) * dq;
  }
}

/**
 * Copy updated temperature values to current array
 * 
 * @param t_current Current temperature array (destination)
 * @param t_new New temperature array (source)
 */
void swap_temperature(double *t_current, double *t_new) {
  for (int i = 0; i < N; i++) {
    t_current[i] = t_new[i];
  }
}

/**
 * Save final results to a CSV file
 * 
 * @param temperature Final temperature array
 * @param dx Cell size
 * @return 0 on success, 1 on failure
 */
int save_final_results(double *temperature, double dx) {
  FILE *file = fopen("conduction-transient.csv", "w");
  if (file == NULL) {
    fprintf(stderr, "Error opening final output file\n");
    return 1;
  }
  
  for (int i = 0; i < N; i++) {
    double x = X0 + (i + 0.5) * dx;  // Cell-center coordinate
    fprintf(file, "%g,%g\n", x, temperature[i]);
  }
  
  fclose(file);
  return 0;
}

/**
 * Run the transient heat conduction simulation
 * 
 * @return 0 on success, non-zero on failure
 */
int run_simulation() {
  // Numerical parameters
  const double dx   = L0/N;          // Cell size
  const double eps  = 0.1;           // Half-width of the initial "rectangle"
  
  // Time-stepping control (explicit)
  // Stability condition for 1D heat eq.: dt <= dx^2/2
  double dt        = 0.5 * dx * dx;  // Time step
  double tmax      = 1.0;            // Final time
  double tprint    = 0.1;            // Print interval for console output
  double tsnap     = 0.1;            // Snapshot interval for file output
  double nextPrint = 0.0;
  double nextSnap  = 0.0;

  // Create output directory
  if (create_output_directory() != 0) {
    return 1;
  }
  
  // Allocate arrays
  double *T  = (double*) calloc(N, sizeof(double)); // Temperature
  double *Tn = (double*) calloc(N, sizeof(double)); // Updated temperature
  double *q  = (double*) calloc(N+1, sizeof(double)); // Flux at interfaces
  
  if (T == NULL || Tn == NULL || q == NULL) {
    fprintf(stderr, "Error allocating memory\n");
    free(T);
    free(Tn);
    free(q);
    return 1;
  }

  // Initialize temperature field
  initialize_temperature(T, dx, eps);

  // Time integration
  double time = 0.0;
  while (time < tmax) {
    // Determine the time step for this iteration
    double current_dt = dt;
    
    // If next time step would cross a snapshot time, adjust dt to hit it exactly
    // !Highlight: this will be done automatically in Basilisk
    if (time < nextSnap && time + dt > nextSnap) {
      current_dt = nextSnap - time;
    }
    
    // Output data to console if print time is reached
    if (time >= nextPrint) {
      print_temperature(T, dx, time);
      nextPrint += tprint;
    }

    // Save snapshot if snapshot time is reached
    if (fabs(time - nextSnap) < 1e-10) {
      save_snapshot(T, dx, time);
      nextSnap += tsnap;
    }

    // Compute fluxes at cell interfaces
    compute_fluxes(T, q, dx);

    // Update temperature field
    update_temperature(T, Tn, q, dx, current_dt);

    // Swap arrays for next time step
    swap_temperature(T, Tn);

    // Advance time
    time += current_dt;
  }

  // Save final results
  int status = save_final_results(T, dx);

  // Free memory
  free(T);
  free(Tn);
  free(q);

  return status;
}

int main() {
  int status = run_simulation();
  return status;
}