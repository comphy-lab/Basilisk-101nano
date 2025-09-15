#include "axi.h"                  // Axisymmetric formulation (r,z) coordinates
#include "navier-stokes/centered.h"  // Main Navier-Stokes solver

/**
 * # Taylor-Culick Viscoelastic Flow Simulation
 *
 * This simulation models the capillary-driven retraction of a fluid filament,
 * known as Taylor-Culick flow, for viscoelastic fluids. 
 * 
 * ## Physical Background
 * 
 * When a thin cylindrical fluid filament with free surface is suddenly released,
 * surface tension causes it to retract. For Newtonian fluids, the retraction velocity
 * reaches a constant value (Taylor-Culick velocity). For viscoelastic fluids, the 
 * retraction dynamics are more complex due to elastic effects.
 * 
 * ## Viscoelastic Fluid Behavior
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
 * ## Exercise 1:
 * After running the simulation, identify different phases of the retraction process.
 * How does the retraction velocity change over time compared to what you'd expect
 * for a Newtonian fluid?
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
 * ### Key Non-Dimensional Numbers:
 * - **Ohnesorge number (Oh)**: Relates viscous forces to inertial and surface tension forces
 * - **Deborah number (De)**: Ratio of material relaxation time to flow characteristic time
 * - **Elasticity number (Ec)**: Measures elastic contribution relative to viscous effects
 * 
 * ### Key Numerical Parameters:
 * - **FILTERED**: Enable density and viscosity jump smoothing for numerical stability
 * - **tsnap**: Time interval between snapshots (default: 1e0)
 * - **fErr**: Error tolerance for volume fraction tracking (1e-3)
 * - **KErr**: Error tolerance for curvature calculation (1e-6)
 * - **VelErr**: Error tolerance for velocity field (1e-3)
 * - **AErr**: Error tolerance for conformation tensor (1e-3)
 * - **Ldomain**: Domain size in characteristic lengths (40.0)
 * - **MAXlevel**: Maximum grid refinement level (10)
 * 
 * ### Exercise 2:
 * Try changing fErr to see how it affects the interface resolution and
 * computational cost. Smaller values give sharper interfaces but require more cells.
 * 
 * ### Exercise 3:
 * Experiment with different Deborah numbers (De) to see how elasticity affects
 * the retraction dynamics. Try De = 0 (Newtonian), De = 1 (moderate elasticity),
 * and De = 10 (high elasticity).
 */
#define FILTERED 1                         // Smear density and viscosity jumps for stability
#include "two-phaseVE.h"                   // Two-phase flow with viscoelasticity
#include "navier-stokes/conserving.h"      // Conservative form of N-S for better mass conservation
#include "tension.h"                       // Surface tension implementation

#define tsnap (1e0)                       // Time interval between output snapshots
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                        // Error tolerance in volume fraction f
#define KErr (1e-6)                        // Error tolerance in curvature calculation
#define VelErr (1e-3)                      // Error tolerances in velocity field
#define AErr (1e-3)                        // Error tolerances in conformation tensor


#define hole0 (1e0)                        // Initial hole size
#define h0 (1e0)                           // Initial filament thickness

/**
 * ## Global Parameters
 * 
 * Here we define the key dimensionless parameters that control the physics:
 * - De: Deborah number (elasticity)
 * - Ohf: Ohnesorge number for fluid 1 (viscosity effect)
 * - Ohs: Ohnesorge number for fluid 2 (surrounding medium)
 * - Ec: Elasticity parameter (polymer concentration)
 * - RhoR: Density ratio between fluids
 */
double De, Ohf, Ohs, Ec, RhoR, tmax, Ldomain;
int MAXlevel;
char nameOut[256], dumpFile[] = "dump"; // Output file name variables

int main(int argc, char const *argv[]) {
  dtmax = 1e-5;                 // Maximum allowed timestep
  Ohs = 1e-5;                   // Ohnesorge number of surrounding fluid
  RhoR = 1e-3;                  // Density ratio (surrounding/filament)

  Ohf = 0.05;                   // Ohnesorge number of filament fluid
  
  Ec = 1e-2;                    // Elasticity parameter
  De = 1e30;                    // Deborah number (relaxation time)
  MAXlevel = 10;                // Maximum grid refinement level

  tmax = 20.0;                  // Maximum simulation time
  Ldomain = 40.0;               // Domain size

  L0=Ldomain;                   // Domain length
  X0=0.0; Y0=0.;                // Domain origin
  init_grid (1 << (10));        // Initialize grid at level 10

  // Create output directory
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  // Set fluid properties
  rho1 = 1.000; mu1 = Ohf;      // Filament fluid properties
  rho2 = RhoR; mu2 = Ohs;       // Surrounding fluid properties

  // Set viscoelastic properties
  G1 = Ec; G2 = 0.0;            // Elastic modulus (G1 for filament)
  lambda1 = De; lambda2 = 0.0;  // Relaxation time (lambda1 for filament)

  f.sigma = 1.0;                // Surface tension coefficient

  // Print simulation parameters
  fprintf(ferr, "Capillary Scaling: Level %d tmax %g. Ohf %3.2e, Ec %3.2e, De %3.2f\n", MAXlevel, tmax, Ohf, Ec, De);
  run();
}

/**
 * ## Initial Condition
 * 
 * This event sets up the initial configuration of the simulation.
 * 
 * We initialize a cylindrical filament with radius h0/2 ending at a 
 * hemispherical cap. The interface is defined using the distance function
 * in the fraction() function, which sets the volume fraction field f.
 * 
 * ### Understanding the Initial Geometry:
 * - For y < hole0+(h0/2.0): We define a hemisphere
 * - Otherwise: We define a cylinder with radius h0/2
 * 
 * ### Exercise 4:
 * Try modifying the initial geometry by changing hole0 and h0.
 * How does the initial shape affect the retraction dynamics?
 */
event init(t = 0){
  if(!restore (file = dumpFile)){
    // Refine grid near the interface
    refine(x<(h0/2.0)+0.025 && y < hole0+(h0/2.0)+0.025 && level<MAXlevel);
    
    // Initialize the interface geometry (cylindrical filament with hemispherical cap)
    fraction(f, y < hole0+(h0/2.0) ? sq(h0/2.0)-(sq(x)+sq(y-(h0/2.0)-hole0)) : (h0/2.0)-x);
  }
  // mask (x > 10.0 ? right : none);
}

/**
 * ## Adaptive Mesh Refinement
 * 
 * This event runs after each timestep to adapt the mesh based on solution features.
 * Adaptation is crucial for efficiency - we want fine cells only where needed.
 * 
 * ### How Adaptive Refinement Works:
 * 1. Basilisk analyzes solution gradients in each cell
 * 2. Cells with steep gradients are divided (refined)
 * 3. Cells with small gradients are merged (coarsened)
 * 4. This balances accuracy and computational cost
 * 
 * ### Refinement Criteria:
 * 1. Interface location (volume fraction gradient)
 * 2. Velocity field gradients 
 * 3. Conformation tensor gradients
 * 4. Interface curvature
 * 
 * ### Exercise 5:
 * Modify the error thresholds (fErr, VelErr, etc.) and observe how
 * they affect the mesh refinement and solution accuracy. What's the
 * trade-off between accuracy and computational cost?
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
 * This event creates simulation output at regular intervals for:
 * 1. Restart capabilities (to continue interrupted simulations)
 * 2. Visualization and analysis
 * 
 * ### Output Files:
 * - **dump**: Contains complete simulation state for restarts
 * - **intermediate/snapshot-X.XXXX**: Contains simulation state at time X.XXXX
 * 
 * ### Visualization Methods:
 * You can use Basilisk's built-in visualization tools or external tools:
 * - Basilisk View: For quick visualization (`bview`)
 * - ParaView: For detailed 3D visualization
 * - Python with matplotlib: For custom analysis
 * 
 * ### Exercise 6:
 * Try changing tsnap to capture more or fewer frames.
 * Smaller tsnap values give smoother animations but larger data size.
 */
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);                  // Create/update restart file
  snprintf (nameOut, sizeof(nameOut), "intermediate/snapshot-%5.4f", t);  // Create snapshot files
  dump(file=nameOut);
}

/**
 * ## Simulation Completion
 * 
 * This event runs at the end of the simulation (t = end) and prints
 * a summary of the simulation parameters to standard error.
 * 
 * ### Important Dimensionless Numbers:
 * - **Wi** (Weissenberg number): Measures elastic effects in the flow
 * - **El** (Elasticity number): Ratio of elastic to inertial forces
 * - **Oh** (Ohnesorge number): Ratio of viscous to inertial and capillary forces
 * - **Bo** (Bond number): Ratio of gravitational to surface tension forces
 * 
 * ### Exercise 7:
 * After completing a simulation, calculate the Taylor-Culick velocity
 * (v = sqrt(2*sigma/(rho*h))) and compare with your simulation results.
 * How do viscoelastic effects modify this theoretical prediction?
 */
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, De %2.1e, Ec %2.1e, Ohf %2.1e, Ohs %2.1e\n", 
            MAXlevel, De, Ec, Ohf, Ohs);
}

/**
 * ## Simulation Logging and Monitoring
 * 
 * This event runs after each timestep to:
 * 1. Calculate the total kinetic energy in the system
 * 2. Log simulation progress and diagnostics
 * 3. Check for simulation stability and stop if necessary
 * 
 * ### Kinetic Energy Calculation:
 * In axisymmetric simulations, volume integrals must include the 2πy factor
 * to account for the cylindrical geometry, where y is the distance from the axis.
 * 
 * ### Stability Monitoring:
 * The simulation automatically detects and stops for:
 * - Energy blow-up (indicates numerical instability)
 * - Energy decay to near-zero (indicates flow has stopped)
 * 
 * ### Exercise 8:
 * Add code to calculate and log additional quantities such as:
 * - Maximum velocity in the domain
 * - Filament length over time
 * - Elastic energy stored in the polymers
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
      fprintf(ferr, "Level %d, De %2.1e, Ec %2.1e, Ohf %2.1e, Ohs %2.1e\n", 
            MAXlevel, De, Ec, Ohf, Ohs);
      fprintf (ferr, "Wi El Oh i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, De %2.1e, Ec %2.1e, Ohf %2.1e, Ohs %2.1e\n", 
            MAXlevel, De, Ec, Ohf, Ohs);
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
     * 
     * ### Exercise 9:
     * Modify the stability criteria thresholds. What happens if you set
     * them too strictly or too loosely?
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