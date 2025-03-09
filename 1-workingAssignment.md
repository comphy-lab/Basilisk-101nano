# Heat Conduction Simulation: From Vanilla C to Basilisk C

## Introduction

This working class will guide you through implementing heat conduction simulations using both vanilla C and the Basilisk framework. We'll progress from simple steady-state problems to more complex transient simulations in multiple dimensions, highlighting the advantages of Basilisk for computational fluid dynamics (CFD) and heat transfer applications.

> [!info] Learning Objectives
> - Implement basic heat conduction solvers in vanilla C
> - Transition to equivalent implementations in Basilisk C
> - Understand the advantages of Basilisk's domain-specific language
> - Develop simulations with increasing complexity (1D → 2D → axisymmetric)
> - Appreciate how Basilisk simplifies numerical methods implementation

## Part 1: Heat Conduction in Vanilla C

### 1.1 Steady-State Heat Conduction (0-conduction-simple.c)

We'll begin with a fundamental 1D steady-state heat conduction problem in vanilla C.

> [!tldr] Problem Statement
> Solve the steady-state heat equation:
> $$\frac{d^2 T}{dx^2} = 0$$
> Subject to Dirichlet boundary conditions:
> - $T(0) = 0$
> - $T(1) = 1$
>
> The analytical solution is the linear profile $T(x) = x$.

Let's examine the key elements of this implementation:

```c
// Domain discretization
#define N        11         // Number of grid points
#define MAX_ITER 10000      // Maximum number of iterations
#define TOL      1e-10      // Convergence tolerance

// Core iterative solver using Gauss-Seidel method
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
```

This implementation uses the Gauss-Seidel iterative method with a central difference approximation for the second derivative:

$$T_{i} \approx \frac{T_{i-1} + T_{i+1}}{2}$$

> [!note] Key Implementation Details
> - Arrays for storing current and new temperature values
> - Explicit boundary condition application
> - Manual iteration until convergence
> - Manual memory management
> - Explicit file I/O for results

### 1.2 Transient Heat Conduction (0-conduction-transient.c)

Next, we'll look at a 1D transient heat conduction problem in vanilla C.

> [!tldr] Problem Statement
> Solve the transient heat equation:
> $$\frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2}$$
> With no-flux boundary conditions and an initial condition approximating a Dirac delta function.
>
> The analytical solution is a Gaussian profile that spreads over time:
> $$T(x,t) = \frac{1}{2\sqrt{\pi t}}e^{-x^2/4t}$$

Key aspects of this implementation:

```c
// Compute heat fluxes at cell interfaces
void compute_fluxes(double *temperature, double *flux, double dx) {
  flux[0] = 0.0; // Left boundary (no flux)
  flux[N] = 0.0; // Right boundary (no flux)
  
  for (int i = 0; i < N - 1; i++) {
    flux[i+1] = - (temperature[i+1] - temperature[i]) / dx;
  }
}

// Update temperature field based on fluxes
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
```

> [!note] Key Implementation Details
> - Explicit time-stepping scheme (forward Euler)
> - Manual computation of fluxes at cell interfaces
> - Explicit handling of boundary conditions
> - Stability condition for time step ($dt \leq dx^2/2$)
> - Manual memory management and array allocation
> - Explicit file I/O for snapshots and results

### Exercise 1: Modify the Vanilla C Code

> [!task] Suggested Tasks
> 1. Modify `0-conduction-simple.c` to 2D.
> 2. Modify `0-conduction-transient.c` to use a different initial condition (e.g., a Gaussian profile)


## Part 2: Introduction to Basilisk C

### 2.1 Steady-State Heat Conduction in Basilisk (0-conduction-simple-basilisk.c)

Now let's examine the equivalent steady-state problem implemented in Basilisk C.

```c
#include "grid/cartesian1D.h"
#include "run.h"

// Declare scalar field for temperature
scalar T[];
T[left] = dirichlet(0.0);
T[right] = dirichlet(1.0);

int main() {
  // Domain setup
  L0 = 1.0;     // Domain length
  X0 = 0.0;     // Left boundary
  N = 200;      // Number of cells
  
  // Set the timestep based on stability criterion
  DT = (L0/N)*(L0/N)/2;
  
  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  
  // Run simulation
  run();
}

// Time integration using explicit finite volume method
event marching (i++) {
  foreach() {
    // Proper heat equation time stepping
    T[] += DT*(T[1] - 2*T[] + T[-1])/(Delta*Delta);
  }
}
```

> [!important] Basilisk Advantages
> - Simplified domain and grid setup
> - Declarative boundary conditions
> - Event-based programming model
> - Intuitive stencil operations with `T[1]` and `T[-1]`
> - Automatic handling of grid traversal with `foreach()`
> - Significantly reduced code length with the same functionality

### 2.2 Understanding Basilisk's Event System

Basilisk uses an event-based programming model that separates concerns:

- `event init (t = 0)`: Initialization at t=0
- `event marching/integration (i++)`: Time integration at each step
- `event end (t = end)`: Final operations at the end of the simulation

Events can also be scheduled at specific times or intervals:

```c
// Execute every 0.1 time units
event writingFiles (t = 0.0; t += 0.1; t < 1.0) {
  // Save output
}
```

### Exercise 2: Understand Basilisk Code Structure

> [!task] Tasks
> 1. Compare the structure and length of `0-conduction-simple.c` and `0-conduction-simple-basilisk.c`
> 2. Identify the advantages of the Basilisk implementation
> 3. Modify `0-conduction-simple-basilisk.c` to output the error relative (not absolute) to the analytical solution at each iteration

## Part 3: Transient Heat Conduction in Basilisk

### Exercise 3 Fill-in-the-Blanks Exercise (0-conduction-transient-basilisk.c)

Let's now implement the transient heat equation in Basilisk. Here's a template with some parts removed for you to fill in:

```c
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
```

> [!hint] Hint for K value 
> What is the stability condition for the 1D heat equation with explicit Euler time-stepping? The maximum stable timestep relates to the spatial discretization (Δx) by what factor?

> [!hint] Hint for flux computation (step 2) 
> Consider how to compute the gradient of temperature between adjacent cells. In Basilisk, how would you access the current cell and its left neighbor? Remember that the heat flux is proportional to the negative temperature gradient.

> [!hint] Hint for temperature change rate (step 3) 
> The change in temperature is related to the divergence of the heat flux. How would you compute the difference between the flux entering the cell (from the left) and the flux leaving the cell (to the right)? Which neighboring flux values do you need to use?

> [!hint] Hint for temperature update (step 4) 
> For an explicit Euler time step, how do you update the current value based on the rate of change and the timestep? Remember the basic formula: new_value = old_value + timestep × rate_of_change.

## 4. Enhanced Transient Solution with Diffusion Module

Now let's look at an enhanced implementation using Basilisk's built-in diffusion module.

```c
#include "grid/multigrid1D.h"
#include "run.h"
#include "diffusion.h"

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
  
  // We can use a larger timestep with the implicit solver
  DT = 0.01;
  
  // Boundary conditions
  T[left] = neumann(0.);
  T[right] = neumann(0.);
  
  // Run simulation
  run();
}

// Initialize temperature field
event init (t = 0) {
  foreach()
    T[] = (fabs(x) < EPS) ? 1.0/EPS/2.0 : 0.0; 
}

// Time integration using implicit diffusion solver
event integration (i++) {
  // Get timestep for this iteration
  double dt = dtnext(DT);
  
  // Use the diffusion() function to solve the equation
  face vector D[];
  foreach_face()
    D.x[] = 1.0; // Constant diffusion coefficient of 1.0
  
  diffusion(T, dt, D);
}
```

> [!important] Advantages of Using diffusion.h
> - Implicit time stepping allows larger timesteps
> - Automatically handles complex boundary conditions
> - More stable and efficient for diffusion problems
> - Uses multigrid methods for faster convergence

## Part 5: Advanced Heat Conduction Problems

### 5.1 2D Heat Conduction (0-conduction-2D.c)

Now, let's extend our simulations to two dimensions.

> [!tldr] Problem Statement
> Solve the 2D transient heat equation:
> $$\frac{\partial T}{\partial t} = \nabla^2 T$$
> With Dirichlet boundary conditions on all sides.

Key aspects of the 2D implementation:

```c
#include "run.h"
#include "diffusion.h"

// Declare scalar field for temperature
scalar T[];

int main() {
  // Domain setup
  L0 = 10.0;     // Domain length
  X0 = -L0/2;    // Left boundary
  init_grid (1 << 7);      // Number of cells
  
  // Larger timestep with implicit solver
  DT = 0.01;
  
  // Boundary conditions
  T[top] = dirichlet(1.);
  T[bottom] = dirichlet(0.);
  T[left] = dirichlet(0.);
  T[right] = dirichlet(0.);
  
  // Run simulation
  run();
}

// Time integration using implicit diffusion solver
event integration (i++) {
  double dt = dtnext(DT);
  
  face vector D[];
  foreach_face()
    D.x[] = 1.0; // Constant diffusion coefficient
  
  mgd = diffusion(T, dt, D);
}
```

> [!note] 2D Implementation Notes
> - Similar structure to 1D, but now using 2D grid
> - Uses `init_grid(1 << 7)` for a 128×128 grid
> - Boundary conditions for all four sides
> - Implicit diffusion solver works in any dimension

### 5.2 Heat Conduction in an Annulus (0-conduction-2D-annulus.c)

Let's look at a more complex geometry: heat conduction in an annular domain.

> [!tldr] Annulus Problem
> Solve the transient heat equation in an annular domain:
> $$\frac{\partial T}{\partial t} = \nabla^2 T$$
> With Dirichlet boundary conditions:
> - $T = 1$ on inner boundary (r = 1)
> - $T = 0$ on outer boundary (r = 4)

```c
#include "run.h"
#include "embed.h"
#include "diffusion.h"

scalar T[];

// Define inner and outer radii
#define INNER_RADIUS 1.0
#define OUTER_RADIUS 4.0

event init (t = 0) {
  // Define annular geometry using embedded boundaries
  solid (cs, fs, difference (sq(OUTER_RADIUS) - sq(x) - sq(y),
                           sq(INNER_RADIUS) - sq(x) - sq(y)));
    
  // Set boundary conditions
  T[embed] = dirichlet (sq(x) + sq(y) < sq(INNER_RADIUS + 0.1) ? 1.0 : 0.0);
}
```

> [!important] Embedded Boundaries
> - `embed.h` provides tools for complex geometries
> - `solid()` defines the embedded boundary
> - `difference()` creates an annular region
> - `T[embed]` sets boundary conditions on the embedded boundary

### 5.3 Axisymmetric Heat Conduction (0-conduction-Axi.c)

Finally, let's look at an axisymmetric problem.

> [!tldr] Axisymmetric Problem
> Solve the heat equation in cylindrical coordinates with axisymmetry:
> $$\frac{\partial T}{\partial t} = \frac{1}{r}\frac{\partial}{\partial r}(r\frac{\partial T}{\partial r}) + \frac{\partial^2 T}{\partial z^2}$$

```c
#include "axi.h"
#include "run.h"
#include "diffusion.h"

scalar T[];

int main() {
  // [...setup code...]
  
  // Boundary conditions
  T[top] = dirichlet(0.0);
  T[left] = neumann(exp(-y*y)); // Gaussian flux profile
  T[right] = dirichlet(0.);
  
  run();
}

event integration (i++) {
  double dt = dtnext(DT);
  
  face vector D[];
  foreach_face()
    D.x[] = fm.x[]; // Use metric terms for axisymmetric problems
  
  mgd = diffusion(T, dt, D);
}
```

> [!note] Axisymmetric Features
> - `axi.h` handles the axisymmetric formulation
> - `y` coordinate represents radial direction (r)
> - `fm.x[]` includes metric terms for the axisymmetric Laplacian
> - Neumann boundary at r=0 (left boundary) with Gaussian profile

### Exercise 4: Advanced Heat Conduction Problems

> [!task] Tasks
> 1. Implement the 2D heat conduction problem
> 2. Compare your results with expected temperature profiles
> 3. Modify the annulus problem to use a time-dependent boundary condition
> 4. Explore the axisymmetric problem with different heat sources

## Conclusion

Throughout this exercise, we've progressed from simple 1D steady-state heat conduction problems in vanilla C to complex geometries and dimensions using Basilisk C. The key takeaways include:

1. Basilisk significantly reduces the complexity and verbosity of numerical simulation code
2. The event-based programming model provides a clean separation of concerns
3. Basilisk's domain-specific language makes numerical methods more intuitive
4. Built-in modules like `diffusion.h` and `embed.h` provide powerful, optimized solvers
5. Extending to higher dimensions and complex geometries is straightforward

> [!tip] Next Steps
> - Explore other physical phenomena like advection-diffusion
> - Implement variable diffusion coefficients
> - Combine with Navier-Stokes solvers for convection problems
> - Use visualization tools like Basilisk View to create animations

## References and Resources

- [Basilisk Documentation](http://basilisk.fr/Front%20Page)
- [Embedded Boundaries in Basilisk](http://basilisk.fr/src/embed.h)
- [Diffusion Module Documentation](http://basilisk.fr/src/diffusion.h)

> [!question] Need Help?
> If you encounter issues with any of the exercises, please check:
> - Code syntax and structure
> - Boundary condition implementation
> - Event scheduling and dependencies
> - Numerical stability considerations