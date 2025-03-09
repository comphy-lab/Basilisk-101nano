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

int main() {
    double T[N], T_new[N];  // Temperature arrays (current and next iteration)
    double dx = 1.0 / (N - 1);  // Grid spacing
    double error;           // Maximum error in current iteration
    int i, iter;            // Loop counters

    // Initialize temperature array with zeros
    // Apply boundary conditions: T(0) = 0, T(N-1) = 1
    for (i = 0; i < N; i++) {
        T[i] = 0.0;
    }
    T[0]     = 0.0;  // Left boundary condition
    T[N - 1] = 1.0;  // Right boundary condition

    // Iterative scheme to solve the equation
    for (iter = 0; iter < MAX_ITER; iter++) {
        // Apply boundary conditions to new array
        T_new[0]     = 0.0;
        T_new[N - 1] = 1.0;

        error = 0.0;
        // Update interior points using central difference approximation
        for (i = 1; i < N - 1; i++) {
            T_new[i] = 0.5 * (T[i - 1] + T[i + 1]);  // Central difference scheme
            double diff = fabs(T_new[i] - T[i]);
            // Track maximum error for convergence check
            if (diff > error) {
                error = diff;
            }
        }

        // Copy the updated solution to T for next iteration
        for (i = 0; i < N; i++) {
            T[i] = T_new[i];
        }

        // Check for convergence based on maximum temperature change
        if (error < TOL) {
            break;
        }
    }

    // Save results to a CSV file with two columns: x and T
    FILE *file = fopen("conduction-simple.csv", "w");
    for (i = 0; i < N; i++) {
        double x = i * dx;  // Physical coordinate
        fprintf(file, "%g,%g\n", x, T[i]);
    }
    fclose(file);

    printf("Number of iterations: %d\n", iter);

    return 0;
}