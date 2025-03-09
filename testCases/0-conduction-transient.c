#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

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

int main() {
    // Numerical parameters
    const double L0   = 10.0;         // Domain length (from -L0/2 to +L0/2)
    const double X0   = -L0/2;        // Left boundary
    const int    N    = 200;          // Number of cells
    const double dx   = L0/(N);       // Cell size
    const double eps  = 0.1;          // Half-width of the initial "rectangle"
    
    // Time-stepping control (explicit)
    // Stability condition for 1D heat eq.: dt <= dx^2/2
    double dt        = (L0/N)*(L0/N)/2;  // Same as 0.5 * dx*dx
    double tmax      = 1.0;              // Final time
    double tprint    = 0.1;              // Print interval for console output
    double tsnap     = 0.1;              // Snapshot interval for file output
    double nextPrint = 0.0;
    double nextSnap  = 0.0;

    // Create intermediate directory for snapshot files
    struct stat st = {0};
    if (stat("intermediate", &st) == -1) {
        if (mkdir("intermediate", 0700) == -1) {
            fprintf(stderr, "Error creating intermediate directory: %s\n", strerror(errno));
            return 1;
        }
    }
    
    // Arrays
    double *T  = (double*) calloc(N, sizeof(double)); // Temperature
    double *Tn = (double*) calloc(N, sizeof(double)); // Updated temperature
    double *q  = (double*) calloc(N+1, sizeof(double)); 
    // q[i+0.5]: flux at interface between cells i and i+1

    // Initialize T with a "Dirac delta" approximated by a thin rectangle
    // T[i] = 1/(2*eps) for |x|<eps, else 0
    // Ensures total integral of T dx = 1
    for (int i = 0; i < N; i++) {
        double x = X0 + (i + 0.5)*dx; // Cell-center coordinate
        if (fabs(x) < eps) {
            T[i] = 1.0/(2.0*eps);
        } else {
            T[i] = 0.0;
        }
    }

    // Time integration
    double time = 0.0;
    while (time < tmax) {
        // Determine the time step to use for this iteration
        double current_dt = dt;
        
        // If next time step would cross a snapshot time, adjust dt to hit it exactly
        if (time < nextSnap && time + dt > nextSnap) {
            current_dt = nextSnap - time;
        }
        
        // If next print is reached or exceeded, output data to console
        if (time >= nextPrint) {
            for (int i = 0; i < N; i++) {
                double x = X0 + (i + 0.5)*dx;
                // Print x, T[i], time
                printf("%g  %g  %g\n", x, T[i], time);
            }
            printf("\n\n");
            nextPrint += tprint;
        }

        // If next snapshot is reached or at exact tsnap interval, save data to file
        if (fabs(time - nextSnap) < 1e-10) {
            char filename[100];
            sprintf(filename, "intermediate/snapshot-%5.4f.csv", time);
            FILE *snap_file = fopen(filename, "w");
            if (snap_file == NULL) {
                fprintf(stderr, "Error opening file %s\n", filename);
            } else {
                for (int i = 0; i < N; i++) {
                    double x = X0 + (i + 0.5)*dx;
                    fprintf(snap_file, "%g,%g\n", x, T[i]);
                }
                fclose(snap_file);
                printf("Saved snapshot at time = %5.4f\n", time);
            }
            nextSnap += tsnap;
        }

        // Compute fluxes at cell interfaces using the formula:
        // q[i+0.5] = - (T[i+1] - T[i])/dx
        q[0]   = 0.0; // Left boundary (no flux)
        q[N]   = 0.0; // Right boundary (no flux)
        for (int i = 0; i < N - 1; i++) {
            q[i+1] = - (T[i+1] - T[i]) / dx;
        }

        // Update T using flux divergence:
        // T^{n+1}_i = T^{n}_i - (dt/dx)*(q[i+0.5] - q[i-0.5])
        for (int i = 0; i < N; i++) {
            double dq;
            if (i == 0) {
                // leftmost cell => flux difference is q[0.5] - q[-0.5], but q[-0.5]=0
                dq = q[1];
            } else if (i == N-1) {
                // rightmost cell => flux difference is q[N-0.5] - q[N-1.5],
                // but q[N+0.5]=0 => effectively just the difference inside domain
                dq = - q[N-1];
            } else {
                dq = q[i+1] - q[i];
            }
            Tn[i] = T[i] - (current_dt/dx)*dq;
        }

        // Swap old and new
        for (int i = 0; i < N; i++) {
            T[i] = Tn[i];
        }

        // Increase time
        time += current_dt;
    }

    // Save final results to a CSV file with two columns: x and T
    FILE *file = fopen("conduction-transient.csv", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening final output file\n");
        return 1;
    }
    
    for (int i = 0; i < N; i++) {
        double x = X0 + (i + 0.5)*dx;  // Cell-center coordinate
        fprintf(file, "%g,%g\n", x, T[i]);
    }
    fclose(file);

    // Free memory
    free(T);
    free(Tn);
    free(q);

    return 0;
}