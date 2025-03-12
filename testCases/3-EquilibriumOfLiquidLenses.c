/**
# Equilibrium of Liquid Lenses in a Three-Phase System

This simulation models the equilibrium shape of a liquid lens at the interface 
between two immiscible fluids. The system consists of three phases:
1. Dispersed phase (FC-40, labeled as phase 1)
2. Water drop (labeled as phase 2)
3. Air (labeled as phase 3)

This physical configuration is relevant to many natural and industrial processes 
including emulsions, droplet coalescence, and interface science.

## Physics and Dimensionless Numbers

The simulation uses several key dimensionless numbers:
- **Ohnesorge number (Oh)**: Relates viscous forces to inertial and surface tension forces
  Oh = μ/√(ρσL), used here for all three phases
  - Oh₁ = 0.5 (dispersed phase)
  - Oh₂ = 0.5 (water drop)
  - Oh₃ = 0.005 (air)
- **Bond number (Bo)**: Relates gravitational forces to surface tension forces
  Bo = ρ₁gR²/σ₁₃ ≈ 1.0, where R is the water drop radius
- **Density ratios**: 
  - RHO21 = ρ₂/ρ₁ ≈ 1/1.85 (water/dispersed phase)
  - RHO31 = ρ₃/ρ₁ ≈ 0.001 (air/dispersed phase)

All dimensionless numbers are based on the properties of the dispersed phase (1)
and the radius of the water drop (2).

## Surface Tension
The surface tension between phases is controlled by:
- f1.sigma: Surface tension coefficient between dispersed phase and air (1-3), set to 1.0 (reference value)
- f2.sigma: Surface tension coefficient between dispersed phase and water drop (1-2), set to 52/16 ≈ 3.25

For an ideal precursor film, these values satisfy the physical constraint:
σ₂₃ = σ₁₃ + σ₁₂

Where σ₂₃ is the effective surface tension between water drop and air (2-3).
This relationship is known as the Neumann triangle condition and represents
the assumption of an ideal precursor film.

## Governing Equations

The dimensionless Navier-Stokes equation being solved is:

Dū/Dt = (1/ρ̄)[-∇p + ∇·(2Oh·D̄) + κ̄δₛn̂] + Bo·ĝ

Where:
- ū is the dimensionless velocity
- p is pressure
- D̄ is the rate-of-strain tensor
- κ̄ is the dimensionless curvature
- δₛ is the interface delta function
- n̂ is the unit normal to the interface
- ĝ is the unit gravity vector

## Implementation Details
- Axisymmetric framework (axi.h) for computational efficiency
- Adaptive mesh refinement for accurate tracking of interfaces
- Careful monitoring of kinetic energy for convergence to equilibrium state
- Volume-Of-Fluid (VOF) interface tracking with the three-phase.h header

## Initial Configuration
A water drop is initialized slightly below the interface between the dispersed phase 
and air. Under the influence of gravity and surface tension, the system will evolve 
to minimize energy and reach equilibrium.

To run this simulation: `CFLAGS=-DDISPLAY=-1 make 3-EquilibriumOfLiquidLenses.tst`
*/

// 1 is the dispersed phase, 2 is water drop, and 3 is air

// Core libraries for simulation
#include "axi.h"                                               // Axisymmetric coordinates
#include "navier-stokes/centered.h"                           // Centered finite-volume Navier-Stokes solver
#define FILTERED                                             // Enable property filtering for stability
#include "three-phase.h"                                     // Three-phase interface tracking
#include "tension.h"                                         // Surface tension model

// Grid refinement levels - control accuracy vs computational cost
#define MAXlevel 8                                              // maximum level
#define MINlevel 3                                              // maximum level

// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-6)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define D2cErr (1e-2)                            // error tolerances in velocity

// Physical domain size (dimensionless)
#define Ldomain 6.0                            // Dimension of the domain

// Dimensionless parameters - configurable for different scenarios
double Oh1, Oh2, Oh3, tmax;                    // Ohnesorge numbers and simulation time
double RHO21, RHO31;                           // Density ratios ρ₂/ρ₁ and ρ₃/ρ₁
double Bo;                                     // Bond number (gravity/surface tension)

#define tsnap (0.10)                           // Time interval between snapshots

/**
## Main Function
*/
int main(int argc, char const *argv[]) {
  
  // Set dimensionless parameters based on physical considerations
  Oh1 = 0.5;                            // Ohnesorge number for dispersed phase: μ₁/√(ρ₁σ₁₃R)
  Oh2 = Oh1;                            // Ohnesorge number for water drop: μ₂/√(ρ₁σ₁₃R)
  Oh3 = 0.005;                          // Ohnesorge number for air: μ₃/√(ρ₁σ₁₃R), 100x less viscous

  Bo = 1.0;                             // Bond number based on dispersed phase: ρ₁gR²/σ₁₃

  RHO21 = 1.0/1.85;                     // Density ratio: water/FC-40 ≈ 1/1.85
  RHO31 = 0.001;                        // Density ratio: air/FC-40 ≈ 0.001

  tmax = 100;                           // Maximum simulation time (ensures equilibrium is reached)

  // Setup the computational domain
  L0=Ldomain;                           // Domain size
  X0=-L0/2.; Y0=0.;                     // Domain origin (axisymmetric, so y=0 is axis of symmetry)
  init_grid (1 << (6));                 // Initialize grid with 2^6 = 64 cells in each direction

  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  // Set physical properties for each phase
  // Note: All properties are dimensionless, normalized by appropriate scales
  rho1 = 1.0; mu1 = Oh1;                // Dispersed phase (reference density and viscosity)
  rho2 = RHO21; mu2 = Oh2;              // Water drop properties
  rho3 = RHO31; mu3 = Oh3;              // Air properties

  // Set surface tension coefficients between phase pairs
  f1.sigma = 1.0;                       // Surface tension: dispersed phase-air (σ₁₃, reference value)
  f2.sigma = 52/16;                     // Surface tension: dispersed phase-water drop (σ₁₂ ≈ 3.25)
  
  // Note: The effective surface tension between water and air (σ₂₃) will be:
  // σ₂₃ = σ₁₃ + σ₁₂ = 1.0 + 52/16 = 4.25
  // This satisfies the ideal precursor film assumption (Neumann triangle)

  // Report simulation parameters
  fprintf(ferr, "tmax %g. surface tension coefficient of dispersed phase -- air %3.2e surface tension coefficient of dispersed phase -- water drop  %3.2e\n", tmax, f1.sigma, f2.sigma);

  run();                                // Start the simulation

}

/**
## Gravitational Acceleration
Add gravity as a constant body force in the x-direction (since x is vertical in axisymmetric coordinates).
The Bond number Bo represents the ratio of gravitational to surface tension forces.
*/
event acceleration(i++) {
  face vector av = a;
  foreach_face(x){
    av.x[] -= Bo;                       // Apply gravitational acceleration proportional to Bond number
  }
}

/**
## Initial Conditions
Setup the initial configuration of the three-phase system with a water drop (phase 2)
positioned slightly below the horizontal interface between the dispersed phase (phase 1)
and air (phase 3). This non-equilibrium configuration will evolve under the combined
effects of gravity and surface tension.
*/
event init(t = 0){

  // Try to restore from a previous simulation state, or initialize from scratch
  if(!restore (file = "dump")){
    // Refine mesh around the water drop for better resolution
    refine(sq(x+1.025) + sq(y) < sq(1.2) && level < MAXlevel);
    
    // Initialize the water drop (phase 2) as a circle slightly below the interface
    // Circle equation: (x+1.025)² + y² = 1, with radius=1 and center at (-1.025,0)
    // The offset of 0.025 below the interface promotes interaction with the interface
    fraction(f2, 1.0 - sq(x+1.025) - sq(y));
    
    // Initialize the dispersed phase (phase 1) in the left half of the domain
    // x<0 is dispersed phase, x>0 is air
    fraction(f1, -x);
    
    // Update the boundary conditions
    boundary ((scalar *){f1, f2});
  }

}

/**
## Adaptive Mesh Refinement
Dynamically refine the mesh based on gradients in key fields to improve accuracy
while maintaining computational efficiency. This is essential for accurately resolving
the triple-junction regions where all three phases meet and the surface tension
balance determines the equilibrium contact angles.
*/
event adapt(i++) {

  // Calculate interface curvatures (needed for surface tension and refinement criteria)
  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);                // Curvature of interface 1
  curvature(f2, KAPPA2);                // Curvature of interface 2

  // Refine mesh based on VOF fractions, velocity fields, and curvatures
  // This ensures high resolution at interfaces and regions of high velocity gradients
  adapt_wavelet ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr},
    MAXlevel, MINlevel);

  // Alternative approach: Include viscous dissipation in refinement criteria
  // This section is commented out but can be enabled for more detailed refinement
  // The viscous dissipation function provides additional information about energy
  // dissipation regions that may require higher resolution
  /*
  scalar KAPPA1[], KAPPA2[], D2c[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  foreach(){
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/y);
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    D2c[] = f1[]*sqrt(sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
  }
  boundary ((scalar *){D2c, KAPPA1, KAPPA2});
  adapt_wavelet ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2, D2c},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr, D2cErr},
    MAXlevel, MINlevel);
  */
}

/**
## Output Snapshots
Save simulation state at regular intervals for visualization and analysis.
This allows tracking the evolution of the liquid lens toward its equilibrium shape.
*/
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");                 // Save restart file (latest state)
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);                // Save snapshot at current time
}

/**
## Monitoring and Convergence Check
Log simulation progress and check for convergence to steady state by monitoring
the kinetic energy of the system. As the liquid lens approaches equilibrium,
the kinetic energy will approach zero, indicating that forces are balanced.
*/
event logWriting (i+=10) {
  // Calculate total kinetic energy (KE) of the dispersed phase
  // The kinetic energy is weighted by cell volume (2*pi*y*Delta²) for axisymmetric coordinates
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += f1[]*(sq(u.x[]) + sq(u.y[]))*(2*pi*y)*sq(Delta);
  }
  
  // Write log data to file and screen
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

  /*
  Convergence criterion: Stop when kinetic energy is very small,
  indicating the system has reached equilibrium. At equilibrium, the lens 
  shape satisfies the balance of surface tension forces at the triple junctions
  as described by the Neumann triangle relationship.
  */
  if (i>1e2 && ke < 1e-8){
    fprintf(ferr, "Equilibrium reached at time %g\n", t);
    dump (file = "finalState");
    return 1;                          // Stop simulation when equilibrium is reached
  }
}
