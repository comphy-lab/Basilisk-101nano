/**
 * # 1D Steady Heat Conduction Solver
 * 
 * This program solves the steady-state heat conduction equation in one dimension:
 * 
 * \frac{d^2 T}{dx^2} = 0
 * 
 * Subject to Dirichlet boundary conditions:
 * - T(0) = 0
 * - T(1) = 1
 * 
 * The domain [0, 1] is discretized into N equally spaced points, and the
 * equation is solved using an iterative Gauss-Seidel method with a central
 * difference approximation for the second derivative:
 * 
 * T_{i} \approx \frac{T_{i-1} + T_{i+1}}{2}
 * 
 * The analytical solution is the linear profile T(x) = x.
 */

#include <stdio.h>
#include <math.h>

#define N        11         // Number of grid points
#define MAX_ITER 10000      // Maximum number of iterations
#define TOL      1e-10      // Convergence tolerance

/**
 * Initialize temperature array and apply boundary conditions
 * 
 * @param temperature The temperature array to initialize
 */
void initialize_temperature(double temperature[N]) {
  int i;
  
  // Initialize with zeros
  for (i = 0; i < N; i++) {
    temperature[i] = 0.0;
  }
  
  // Apply boundary conditions
  temperature[0] = 0.0;     // Left boundary condition
  temperature[N - 1] = 1.0; // Right boundary condition
}

/**
 * Apply boundary conditions to temperature array
 * 
 * @param temperature The temperature array to update
 */
void apply_boundary_conditions(double temperature[N]) {
  temperature[0] = 0.0;     // Left boundary condition
  temperature[N - 1] = 1.0; // Right boundary condition
}

/**
 * Update interior points using central difference scheme
 * 
 * @param t_current Current temperature array
 * @param t_new New temperature array to be computed
 * @return Maximum error between current and new values
 */
double update_interior_points(double t_current[N], double t_new[N]) {
  int i;
  double error = 0.0;
  double diff;
  
  for (i = 1; i < N - 1; i++) {
    // Central difference scheme
    t_new[i] = 0.5 * (t_current[i - 1] + t_current[i + 1]);
    
    // Track maximum error for convergence check
    diff = fabs(t_new[i] - t_current[i]);
    if (diff > error) {
      error = diff;
    }
  }
  
  return error;
}

/**
 * Copy new temperature values to current array
 * 
 * @param t_current Current temperature array (destination)
 * @param t_new New temperature array (source)
 */
void update_solution(double t_current[N], double t_new[N]) {
  int i;
  
  for (i = 0; i < N; i++) {
    t_current[i] = t_new[i];
  }
}

/**
 * Write temperature results to a CSV file
 * 
 * @param temperature The final temperature array
 * @param dx Grid spacing
 */
void write_results(double temperature[N], double dx) {
  int i;
  FILE *file = fopen("conduction-simple.csv", "w");
  
  if (file == NULL) {
    fprintf(stderr, "Error opening output file\n");
    return;
  }
  
  for (i = 0; i < N; i++) {
    double x = i * dx;  // Physical coordinate
    fprintf(file, "%g,%g\n", x, temperature[i]);
  }
  
  fclose(file);
}

/**
 * Solve 1D steady heat conduction problem
 * 
 * @return Number of iterations performed
 */
int solve_heat_conduction() {
  double t_current[N];  // Current temperature array
  double t_new[N];      // Next iteration temperature array
  double dx = 1.0 / (N - 1);  // Grid spacing
  double error;         // Maximum error in current iteration
  int iter;             // Loop counter
  
  // Initialize temperature array
  initialize_temperature(t_current);
  
  // Iterative scheme to solve the equation
  for (iter = 0; iter < MAX_ITER; iter++) {
    // Apply boundary conditions to new array
    apply_boundary_conditions(t_new);
    
    // Update interior points using central difference approximation
    error = update_interior_points(t_current, t_new);
    
    // Copy the updated solution
    update_solution(t_current, t_new);
    
    // Check for convergence
    if (error < TOL) {
      break;
    }
  }
  
  // Write results to file
  write_results(t_current, dx);
  
  return iter;
}

int main() {
  int iterations = solve_heat_conduction();
  printf("Number of iterations: %d\n", iterations);
  
  return 0;
}