/**
 * # Soft Impact Dynamics with Viscoelastic Fluids
 *
 * This simulation explores the fascinating physics of soft impact phenomena 
 * using viscoelastic fluids. Applications include droplet impact dynamics, 
 * inkjet printing, coating processes, and soft material interactions.
 *
 * ## Learning Objectives
 *
 * By working with this simulation, you will:
 * 1. Understand how viscoelasticity affects impact dynamics
 * 2. Learn to implement axisymmetric simulations in Basilisk
 * 3. Practice working with dimensionless numbers in fluid mechanics
 * 4. Explore adaptive mesh refinement (AMR) for interface tracking
 * 5. Analyze stability considerations in complex multiphase flows
 *
 * ## Hands-on Exercises (suggested)
 *
 * - Modify the Weber number to observe surface tension effects
 * - Change the Weissenberg number to see how elasticity affects jet formation
 * - Experiment with different initial conditions
 * - Compare scalar vs. tensor viscoelastic models
 */

#include "axi.h"                  // Axisymmetric formulation (r,z) coordinates
#include "navier-stokes/centered.h"  // Main Navier-Stokes solver

/**
 * ## Viscoelastic Model Implementation
 * 
 * Viscoelastic fluids (like polymer solutions) show both viscous and elastic behaviors:
 * - Viscous: Resist flow like honey (dissipates energy)
 * - Elastic: Spring-like memory effects (stores energy)
 *
 * We use the log-conformation approach for numerical stability when simulating viscoelastic fluids.
 * This method transforms the conformation tensor logarithmically to preserve positive-definiteness
 * during the numerical solution, avoiding the "high Weissenberg number problem" that plagues
 * traditional approaches.
 * 
 * Reference: V. Sanjay, Zenodo, DOI: 10.5281/zenodo.14210635 (2024)
 */
// #define _SCALAR // uncomment to use the scalar version of the viscoelastic code
#if !_SCALAR
#include "log-conform-viscoelastic.h"  // Tensor implementation (full anisotropic stresses)
#else 
#include "log-conform-viscoelastic-scalar-2D.h"  // Scalar implementation (faster but less accurate)
#endif

/**
 * ## Simulation Parameters
 * 
 * These parameters control the numerical accuracy and computational demands.
 * Adjust them based on your computational resources and required accuracy.
 * 
 * ### Key Numerical Parameters:
 * - **FILTERED**: Enable density and viscosity jump smoothing for numerical stability
 * - **tsnap**: Time interval between snapshots (default: 1e-2)
 * - **fErr**: Error tolerance for volume fraction tracking (1e-3)
 * - **KErr**: Error tolerance for curvature calculation (1e-6)
 * - **VelErr**: Error tolerance for velocity field (1e-3)
 * - **AErr**: Error tolerance for conformation tensor (1e-3)
 * - **Ldomain**: Domain size in characteristic lengths (4)
 * 
 * ### Exercise:
 * Try changing fErr to see how it affects the interface resolution and
 * computational cost. Smaller values give sharper interfaces but require more cells.
 */
#define FILTERED 1                         // Smear density and viscosity jumps for stability
#include "two-phaseVE.h"                   // Two-phase flow with viscoelasticity
#include "navier-stokes/conserving.h"      // Conservative form of N-S for better mass conservation
#include "tension.h"                       // Surface tension implementation

#define tsnap (1e-2)                       // Time interval between output snapshots
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                        // Error tolerance in volume fraction f
#define KErr (1e-6)                        // Error tolerance in curvature calculation
#define VelErr (1e-3)                      // Error tolerances in velocity field
#define AErr (1e-3)                        // Error tolerances in conformation tensor

// Domain size
#define Ldomain 4

/**
 * ## Dimensionless Numbers
 * 
 * Fluid mechanics often uses dimensionless numbers to characterize flow regimes:
 * 
 * - **Weber (We)**: Ratio of inertia to surface tension
 *   We = ρU²L/σ (high We = inertia dominates, low We = surface tension dominates)
 *   
 * - **Ohnesorge (Oh)**: Ratio of viscous forces to inertial and surface tension forces
 *   Oh = μ/√(ρσL) (high Oh = viscosity important)
 *   
 * - **Weissenberg (Wi)**: Ratio of elastic forces to viscous forces
 *   Wi = λγ̇ (relaxation time × strain rate; high Wi = elastic effects dominant)
 *   
 * - **Elastic number (El)**: Ratio of elastic stresses to viscous stresses
 *   El = Wi/Re (measures relative importance of elasticity vs viscosity)
 *   
 * - **Bond (Bo)**: Ratio of gravitational forces to surface tension
 *   Bo = ρgL²/σ (high Bo = gravity dominates, low Bo = surface tension dominates)
 */
#define We 1.0                             // Weber number (surface tension vs inertia)

/**
 * ## Initial Droplet Shape
 * 
 * This function defines the initial shape of the droplet as a perfect circle.
 * Returns the squared distance from the origin (x²+y²), which is used to:
 * 1. Define the interface location (where R2Drop = 1.0)
 * 2. Control mesh refinement near the interface
 * 
 * Exercise: Modify this function to create elliptical droplets or other shapes
 * and observe how the impact dynamics change.
 */
double R2Drop(double x, double y) {
  return sq(x) + sq(y);
}

/**
 * ## Boundary Conditions
 * 
 * In axisymmetric coordinates (r,z):
 * - left boundary = centerline (axis of symmetry, r=0)
 * - right boundary = far-field boundary
 * - bottom boundary = lower z-boundary
 * - top boundary = upper z-boundary
 * 
 * On the centerline (r=0), we enforce:
 * - No radial velocity (u.t[left] = 0)
 * - No fluid crossing boundary (f[left] = 0)
 * 
 * Other boundaries use default conditions (outflow/symmetry)
 */
u.t[left] = dirichlet(0.0);  // No radial velocity on axis of symmetry
f[left] = dirichlet(0.0);    // No fluid crosses the axis
// all other boundaries use symmetry planes by default

// Global parameters
int MAXlevel;                              // Maximum refinement level
double Oh, Oha, Wi, El, Bond, tmax;        // Dimensionless numbers
char nameOut[80], dumpFile[80];            // Output file names

/**
 * ## Main Function
 * 
 * This function sets up the simulation parameters, initializes the domain, 
 * and runs the simulation. The setup follows this sequence:
 * 
 * 1. Set the timestep limit for stability
 * 2. Configure the domain size and origin
 * 3. Define physical parameters (dimensionless numbers)
 * 4. Initialize the grid
 * 5. Set fluid properties based on dimensionless numbers
 * 6. Run the simulation
 */
int main(int argc, char const *argv[]) {
  // Maximum allowed timestep (important for stability)
  dtmax = 1e-5;                            // BEWARE: Critical for numerical stability!

  // Set up domain size and origin
  L0 = Ldomain;                            // Set domain size
  origin (-L0/2., 0.);                     // Center domain at (-L0/2, 0) - needed for axisymmetric flows
  
  /**
   * ### Parameter Setup
   * 
   * In real research, these values would typically be passed from command line.
   * For this educational example, we use fixed representative values.
   * 
   * Physics interpretation of these parameters:
   * - Wi = 1e30: Extremely elastic fluid (approaches purely elastic limit)
   * - El = 0.1: Moderate elasticity compared to viscosity
   * - Oh = 1e-2: Low viscosity fluid (inertia-dominated)
   * - Bond = 5e-1: Moderate gravity effects
   * 
   * EXERCISE: Try changing these parameters and observe the effect on jet dynamics!
   * For example, try:
   * 1. Wi = 0 (Newtonian fluid)
   * 2. Wi = 10 (moderate elasticity)
   * 3. Oh = 1.0 (viscous-dominated)
   */
  MAXlevel = 10;                           // Maximum refinement level (computational cost ~ 4^MAXlevel)
  Wi = 1e30;                               // Weissenberg number (elasticity)
  El = 0.1;                                // Elastic number 
  Oh = 1e-2;                               // Ohnesorge number (solvent viscosity)
  Bond = 5e-1;                             // Bond number (gravity)
  tmax = 1e0;                              // Maximum simulation time

  // To use command line arguments, uncomment these:
  // MAXlevel = atoi(argv[1]);
  // Wi = atof(argv[2]);  // Use a value of 1e30 to simulate the Wi \to \infty limit
  // El = atof(argv[3]);
  // Oh = atof(argv[4]);
  // Bond = atof(argv[5]);
  // tmax = atof(argv[6]);
  
  // Verify command line arguments (commented out for this example)
  // if (argc < 7){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
  //   return 1;
  // }
  
  /**
   * ### Grid Initialization
   * 
   * We start with a coarse base grid (2^5 = 32 cells per dimension)
   * and will use adaptive refinement to add cells where needed.
   * 
   * This approach is much more efficient than using a uniform fine grid.
   * The base grid resolution is intentionally coarse to save memory.
   */
  init_grid (1 << 5);  // 2^5 = 32 cells baseline grid
  
  /**
   * ### Output Setup
   * 
   * Create a directory for intermediate results and set up the restart filename.
   * Basilisk can generate snapshots at regular intervals and restart files
   * for continuing simulations from breakpoints.
   */
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");

/**
 * ### Physical Properties
 * 
 * Here we define the physical properties of both fluids in the simulation.
 * All properties are non-dimensionalized to make the simulation general.
 * 
 * Using non-dimensional values lets us represent entire classes of flows
 * that have the same dimensionless parameters, regardless of scale.
 * 
 * Property overview:
 * - **rho1, rho2**: Density ratio between liquid and gas (1:0.001)
 * - **mu1, mu2**: Viscosity ratio controlled by Ohnesorge numbers
 * - **lambda1, lambda2**: Relaxation times (polymer has memory, gas does not)
 * - **G1, G2**: Elastic moduli (polymer is elastic, gas is not)
 * - **f.sigma**: Surface tension coefficient (normalized to 1.0)
 * 
 * NOTE: For viscoelastic fluids, the total stress tensor has three contributions:
 * σ = -pI + μs(∇u + ∇u^T) + G(A-I)
 * where:
 * - p is pressure
 * - μs is solvent viscosity
 * - G is elastic modulus
 * - A is the conformation tensor (represents polymer stretching)
 */
  // Density ratio (liquid : gas = 1000:1)
  rho1 = 1., rho2 = 1e-3;                  // Phase 1 is heavy fluid, Phase 2 is light fluid
  
  // Gas phase Ohnesorge number (typically lower than liquid)
  Oha = 2e-2 * Oh;                         // Gas is less viscous than liquid
  
  // Viscosity from Ohnesorge numbers
  mu1 = Oh/sqrt(We), mu2 = Oha/sqrt(We);   // Convert Oh to Reynolds number
  
  // Polymer relaxation time (liquid has elasticity, gas does not)
  lambda1 = Wi; lambda2 = 0.;              // Phase 1 is viscoelastic, Phase 2 is Newtonian
  
  // Elastic modulus (liquid is viscoelastic, gas is not)
  G1 = El; G2 = 0.;                        // Only Phase 1 has elastic stresses

  // Surface tension coefficient (normalized to 1.0)
  f.sigma = 1.0/We;                        // Surface tension scaled by Weber number

  /**
   * ### Numerical Parameters
   * 
   * These parameters control the convergence and stability of the simulation:
   * 
   * - TOLERANCE: Controls precision of the pressure solver
   * - CFL: Courant-Friedrichs-Lewy number (stability criterion)
   *   - CFL < 1 required for explicit methods
   *   - Lower CFL = more stable but slower simulation
   */
  TOLERANCE = 1e-4;                        // Tolerance for iterative pressure solver
  CFL = 1e-1;                              // CFL number (stability criterion)

  // Start the simulation
  run();
}

/**
 * ## Initialization Event
 * 
 * This event runs at t=0 and sets up the initial condition.
 * It either:
 * 1. Restores a previous simulation state from a restart file, or
 * 2. Creates a fresh simulation with a spherical droplet
 * 
 * For the fresh start, we:
 * 1. Refine the mesh near the droplet interface
 * 2. Initialize the volume fraction (f) to create the droplet
 * 3. Set the initial velocity field (-1.0 in the droplet, 0 elsewhere)
 * 
 * The negative x-velocity (-1.0) means the droplet is moving left,
 * toward the axis of symmetry, creating an impact scenario.
 */
event init (t = 0) {
  if (!restore (file = dumpFile)){
    // Refine mesh around interface (where R2Drop ≈ 1)
    refine((R2Drop(x, y) < 1.05) && (level < MAXlevel));
    
    // Initialize volume fraction (f=1 inside droplet, f=0 outside)
    fraction(f, 1. - R2Drop(x, y));
    
    // Set initial velocity field
    foreach() {
      u.x[] = -1.0 * f[];  // Droplet moves left with velocity -1.0
      u.y[] = 0.0;         // No vertical motion initially
    }
  }
}

/**
 * ## Adaptive Mesh Refinement
 * 
 * This event runs after each timestep to adapt the mesh based on solution features.
 * Adaptation is crucial for efficiency - we want fine cells only where needed.
 * 
 * The refinement criteria include:
 * 1. Interface location (volume fraction gradient)
 * 2. Velocity field gradients
 * 3. Conformation tensor gradients
 * 4. Interface curvature
 * 
 * Cells that meet refinement criteria are split, while cells in smooth
 * regions are merged to save computational resources.
 */
event adapt(i++){
  // Calculate interface curvature for refinement criterion
  scalar KAPPA[];
  curvature(f, KAPPA);

  // Adapt mesh based on solution gradients
  #if !_SCALAR
   // For tensor viscoelastic model, refine based on all tensor components
   adapt_wavelet ((scalar *){f, u.x, u.y, conform_p.x.x, conform_p.y.y, conform_p.y.x, conform_qq, KAPPA},
      (double[]){fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr},
      MAXlevel, MAXlevel-6);
  #else
   // For scalar viscoelastic model, refine based on scalar components
   adapt_wavelet ((scalar *){f, u.x, u.y, A11, A22, A12, AThTh, KAPPA},
      (double[]){fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr},
      MAXlevel, MAXlevel-6);
  #endif
}

/**
 * ## Output Generation
 * 
 * This event creates simulation output at regular intervals.
 * It saves:
 * 1. Restart files (for continuing simulations)
 * 2. Snapshot files (for visualization and analysis)
 * 
 * Visualization can be done with tools like ParaView, VisIt, or
 * Basilisk's built-in visualization functions.
 * 
 * EXERCISE: Try changing tsnap to capture more or fewer frames
 * of the simulation for visualization.
 */
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);                  // Create/update restart file
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);  // Create snapshot files
  dump(file=nameOut);
}

/**
 * ## Simulation Completion
 * 
 * This event runs at the end of the simulation (t = end) and prints
 * a summary of the simulation parameters to standard error.
 * 
 * This serves as a record of what parameters were used for this run,
 * which is helpful when running multiple simulations for comparison.
 */
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, Wi %2.1e, El %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", 
            MAXlevel, Wi, El, Oh, Oha, Bond);
}

/**
 * ## Simulation Logging
 * 
 * This event runs after each timestep to:
 * 1. Calculate the total kinetic energy in the system
 * 2. Log simulation progress and diagnostics
 * 3. Check for simulation stability and stop if necessary
 * 
 * The kinetic energy is a good indicator of simulation health:
 * - If it grows too large, the simulation is becoming unstable
 * - If it becomes too small, the flow has essentially stopped
 * 
 * IMPORTANT: For axisymmetric simulations, volume integrals must include
 * the 2πy factor in cylindrical coordinates (y is the distance from the axis).
 */
event logWriting (i++) {
  // Calculate total kinetic energy in the domain
  double ke = 0.;
  foreach (reduction(+:ke)){
    // For axisymmetric coordinates, include 2πy factor for volume element
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  
  // Write log information (only on main process in parallel runs)
  if (pid() == 0) {
    static FILE * fp;
    if (i == 0) {
      // Write header information at start
      fprintf(ferr, "Level %d, Wi %2.1e, El %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", 
              MAXlevel, Wi, El, Oh, Oha, Bond);
      fprintf (ferr, "Wi El Oh i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, Wi %2.1e, El %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", 
              MAXlevel, Wi, El, Oh, Oha, Bond);
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      // Append data in subsequent timesteps
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

    /**
     * ### Stability Checks
     * 
     * These checks monitor the simulation health and stop if:
     * 1. Kinetic energy becomes negative (physical impossibility)
     * 2. Kinetic energy grows too large (numerical instability)
     * 3. Kinetic energy becomes too small (simulation essentially completed)
     * 
     * This prevents wasting computational resources on failed or finished runs.
     */
    
    // Check for negative kinetic energy (physical impossibility, numerical error)
    assert(ke > -1e-10);

    // Check for kinetic energy blowup (simulation instability)
    if (ke > 1e2 && i > 1e1){
      if (pid() == 0){
        fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
        fp = fopen ("log", "a");
        fprintf(fp, "The kinetic energy blew up. Stopping simulation\n");
        fclose(fp);
        dump(file=dumpFile);
        return 1;
      }
    }
    assert(ke < 1e2);
    
    // Check if simulation has essentially stopped
    if (ke < 1e-6 && i > 1e1){
      if (pid() == 0){
        fprintf(ferr, "kinetic energy too small now! Stopping!\n");
        dump(file=dumpFile);
        fp = fopen ("log", "a");
        fprintf(fp, "kinetic energy too small now! Stopping!\n");
        fclose(fp);
        return 1;
      }
    }
  }
}

/**
 * ## Data Analysis and Visualization (Post-Processing)
 * 
 * After the simulation is complete, you can use the snapshot files
 * for analysis and visualization. Suggested approaches:
 *  
 * Write custom post-processing scripts for specific metrics
 *    - Jet height vs. time
 *    - Maximum pressure at impact
 *    - Energy dissipation rate
 * 
 */
