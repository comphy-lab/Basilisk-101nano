/**
# Planar Couette flow of Newtonian Fluid

This code simulates a planar Couette flow with a Newtonian fluid.
It uses a constant viscosity model, which is a simplified version of
the generalized non-Newtonian model in 2-Couette-flow-Bingham.c.

## Mathematical Formulations

For Newtonian fluids, the stress-strain rate relationship is linear:
$$ \tau = 2\mu_0 D_{ij} $$

Where $D_{ij}$ is the deformation tensor:
$D_{ij}=(u_{i,j}+u_{j,i})/2$

## Exact solution in the proposed case

For a Newtonian fluid with constant viscosity $\mu_0$ and a constant 
pressure gradient $-\frac{dp}{dx} = 1$, the velocity profile is parabolic:

$$ u(y) = \frac{1}{2\mu_0}y(h-y) $$

Where $h$ is the height of the domain (1.0 in this case).

# Code
*/
#include "navier-stokes/centered.h"

// Global parameters
char file_name[80];
double mu_0 = 1.0;        // Constant viscosity (Newtonian)
int max_iter = 1e4;       // Maximum iterations
#define DT_MAX (1e-3)     // Maximum timestep

int main(int argc, char const *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s filename [mu_0]\n", argv[0]);
    return 1;
  }
  
  sprintf(file_name, "%s", argv[1]);
  
  // Initialize grid and domain
  init_grid(1<<6);
  L0 = 1.0;
  origin(0.0, 0.0);
  DT = DT_MAX;
  stokes = true;
  TOLERANCE = 1e-5;

/**
 For Newtonian fluid, we only need to set the constant viscosity $\mu_0$
*/
  // Parse command line arguments if provided
  if (argc >= 3) {
    mu_0 = atof(argv[2]);
  }

  // Set boundary conditions
  // Right-left boundaries are periodic
  periodic(right);
  
  // Slip at the top
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  uf.n[top] = neumann(0);
  
  // No slip at the bottom
  u.n[bottom] = dirichlet(0);
  uf.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  
  // Pressure conditions are neumann 0.0
  p[top] = neumann(0);
  pf[top] = neumann(0);
  p[bottom] = neumann(0);
  pf[bottom] = neumann(0);

  run();
}

// un is used to search for a stationary solution
scalar un[];

/**
## Initialization event
*/
event init(t = 0) {
  // Set constant Newtonian viscosity
  const scalar mu_const[] = mu_0;
  mu = mu_const;
  
  /**
   Pressure gradient `mdpdx`
   $$-\frac{dp}{dx} = 1 $$
  */
  const face vector mdpdx[] = {1.0, 0.0};
  a = mdpdx;
  
  // Initialize velocity field at rest
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
    un[] = 0;
  }
  
  dump(file = "start");
}

/**
## Monitoring convergence 
We look for a stationary solution by checking changes in velocity field.
*/
event logfile(i += 500; i <= max_iter) {
  double du = change(u.x, un);
  fprintf(ferr, "i = %d: err = %g\n", i, du);
  
  if (i > 0 && du < 1e-6) {
    dump(file = file_name);
    return 1; /* stop */
  }
  
  if (i == max_iter) {
    dump(file = file_name);
  }
}

/**
## Compare with exact solution
*/
event end(t = end) {
  FILE *fp = fopen("profile.dat", "w");
  fprintf(fp, "# y u_numerical u_exact error\n");
  
  // Sample the velocity profile along the x-midpoint
  foreach_col(0) {
    // Calculate exact solution: u(y) = (1/2μ) * y * (h-y)
    double y = y;
    double u_exact = (1.0/(2.0*mu_0)) * y * (L0-y);
    double error = fabs((u.x[] - u_exact)/u_exact);
    
    fprintf(fp, "%g %g %g %g\n", y, u.x[], u_exact, error);
  }
  
  fclose(fp);
  fprintf(ferr, "Profile data written to 'profile.dat'\n");
}
