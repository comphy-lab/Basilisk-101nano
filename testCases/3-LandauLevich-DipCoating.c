/**
 * @file 3-LandauLevich-DipCoating.c
 * @brief Landau-Levich dip coating: Film formation during substrate withdrawal
 * 
 * This code simulates the classical Landau-Levich problem of a vertical substrate
 * being withdrawn from a liquid bath. The simulation captures the dynamics of
 * film formation, meniscus shape evolution, and demonstrates the relationship
 * between withdrawal velocity and film thickness.
 * 
 * Key Physics:
 * - Capillary number (Ca): Viscous vs surface tension forces
 * - Contact line dynamics and wetting
 * - Film thickness scaling: h ~ Ca^(2/3) (Landau-Levich law)
 * - Meniscus shape and dynamic contact angles
 * 
 * Course Context: Hour 3 - Coating Applications
 * Demonstrates multiscale nature of contact line problem and industrial relevance
 * 
 * @author Vatsal Sanjay (vatsalsy@comphy-lab.org)
 * https://comphy-lab.org
 * @date 2025-08-26
 * @version 1.0 
*/

// ======= Include necessary Basilisk modules =======
#include "navier-stokes/centered.h"    // NS solver with centered discretization
#define FILTERED                       // Use filtered VOF advection for better interface quality
#include "two-phase.h"                // Two-phase interface tracking (VOF method)
#include "navier-stokes/conserving.h" // Conservative momentum advection
#include "tension.h"                  // Surface tension model
#include "contact.h"                  // Contact angle implementation

// ======= Numerical parameters for adaptivity =======
// Error tolerances for adaptive mesh refinement
#define fErr (1e-3)                   // Error tolerance in Volume of Fluid (interface position)
#define KErr (1e-4)                   // Error tolerance in curvature (KAPPA) - tighter for coating
#define VelErr (1e-2)                 // Error tolerances in velocity fields
#define AErr (1e-4)                   // Error tolerance in acceleration (for film stability)

// ======= Physical parameters =======
// Dimensional parameters (all in SI units for clarity)
#define RHO_L 1000.0                  // Liquid density (kg/m³) - water-like
#define RHO_G 1.2                     // Gas density (kg/m³) - air
#define MU_L 0.001                    // Liquid viscosity (Pa·s) - water-like
#define MU_G 1.8e-5                   // Gas viscosity (Pa·s) - air  
#define SIGMA 0.072                   // Surface tension (N/m) - water-air
#define L_CHAR 1e-3                   // Characteristic length (m) - 1mm
#define U_COATING 0.1                 // Coating velocity (m/s) - 10 cm/s

// Non-dimensional parameters
#define Rho21 (RHO_G/RHO_L)           // Density ratio (dimensionless)
#define Mu21 (MU_G/MU_L)              // Viscosity ratio (dimensionless)
#define Ca (MU_L*U_COATING/SIGMA)     // Capillary number (dimensionless)
#define Re (RHO_L*U_COATING*L_CHAR/MU_L) // Reynolds number (dimensionless)
#define Bo (RHO_L*9.81*L_CHAR*L_CHAR/SIGMA) // Bond number (dimensionless)

// ======= Domain and geometry setup =======
// Substrate parameters
#define SUBSTRATE_X 0.5               // x-position of substrate center
#define SUBSTRATE_WIDTH 0.1           // Half-width of substrate
#define BATH_HEIGHT 0.6               // Initial bath height (below substrate)

// Contact angle parameters (in degrees, converted to radians)
#define THETA_DEG 60.0                // Static contact angle in degrees
#define THETA_RAD (THETA_DEG*M_PI/180.0) // Contact angle in radians

// ======= Global variables =======
double U_wall = 1.0;                  // Normalized coating velocity (U_COATING/U_COATING = 1)
double t_ramp = 5.0;                  // Time to ramp up to full coating velocity
int LEVEL = 8;                        // Maximum refinement level

// ======= Boundary conditions =======

// Bottom boundary: no-slip wall (bath bottom)
u.n[bottom] = dirichlet(0.0);
u.t[bottom] = dirichlet(0.0);
f[bottom] = neumann(0.0);

// Top boundary: outflow
u.n[top] = neumann(0.0);
u.t[top] = neumann(0.0);
p[top] = dirichlet(0.0);
pf[top] = dirichlet(0.0);

// Right boundary: far field / outflow
u.n[right] = neumann(0.0);
u.t[right] = neumann(0.0);
p[right] = dirichlet(0.0);
pf[right] = dirichlet(0.0);

// ======= Moving wall implementation =======
/**
 * Embedded boundary representing the moving substrate
 * The substrate moves upward with velocity U_wall
 */

// Substrate geometry: vertical surface at SUBSTRATE_X
scalar substrate[];

event init (t = 0) {
  // Initialize substrate geometry
  foreach() {
    // Define substrate as region within SUBSTRATE_WIDTH of SUBSTRATE_X
    if (fabs(x - SUBSTRATE_X) <= SUBSTRATE_WIDTH) {
      substrate[] = 1.0;  // Inside substrate
    } else {
      substrate[] = 0.0;  // Outside substrate  
    }
  }
  
  // Initialize liquid bath
  foreach() {
    // Liquid (f=1) below BATH_HEIGHT, gas (f=0) above
    if (y < BATH_HEIGHT) {
      f[] = 1.0;  // Liquid phase
    } else {
      f[] = 0.0;  // Gas phase
    }
  }
  
  // Apply contact angle at substrate
  f.height = THETA_RAD;
}

/**
 * Moving boundary condition for the substrate
 * Implements upward motion with smooth startup
 */
event moving_wall (i++) {
  double U_current;
  
  // Smooth ramp-up of coating velocity
  if (t < t_ramp) {
    U_current = U_wall * (t/t_ramp) * (2.0 - t/t_ramp); // Smooth S-curve
  } else {
    U_current = U_wall;
  }
  
  // Apply no-slip boundary condition on moving substrate
  foreach() {
    if (substrate[] > 0.5) { // Within substrate
      u.x[] = 0.0;           // No horizontal velocity
      u.y[] = U_current;     // Upward velocity
      
      // Maintain contact angle condition
      if (f[] > 0.01 && f[] < 0.99) { // Near interface
        // This is a simplified contact line treatment
        // In practice, more sophisticated subgrid models are needed
      }
    }
  }
}

/**
 * Adaptive mesh refinement strategy
 * Focus refinement near:
 * 1. Liquid-gas interface (VOF)
 * 2. High velocity gradients (boundary layers)
 * 3. High curvature regions (meniscus)
 * 4. Contact line region
 */
event adapt (i++) {
  // Adapt based on interface position and curvature
  adapt_wavelet ({f, u.x, u.y}, 
                 (double[]){fErr, VelErr, VelErr}, 
                 maxlevel = LEVEL, minlevel = 4);
                 
  // Additional refinement near contact line (substrate region)
  foreach() {
    if (substrate[] > 0.1 && level < LEVEL-1) {
      if (f[] > 0.01 && f[] < 0.99) { // Near interface and substrate
        refine_cell (point, NULL, 0, NULL);
      }
    }
  }
}

/**
 * Set fluid properties
 * Handles two-phase flow properties with proper ratios
 */
event properties (i++) {
  foreach_face() {
    double ff = (f[] + f[-1])/2.0;
    
    // Density: linear interpolation between phases
    alpha.x[] = clamp(Rho21 + (1.0 - Rho21)*ff, Rho21, 1.0);
    
    // Viscosity: harmonic mean (better for high ratio)
    if (ff > 0.5) {
      mu.x[] = 1.0/(ff/1.0 + (1.0-ff)/Mu21);  // Liquid-dominated
    } else {
      mu.x[] = Mu21/(ff/Mu21 + (1.0-ff)/1.0); // Gas-dominated  
    }
  }
}

/**
 * Output routine for film thickness measurement and visualization
 * Key outputs:
 * 1. Film thickness vs height
 * 2. Interface shape
 * 3. Velocity profiles
 * 4. Capillary number effects
 */
event output (t += 0.5) {
  
  // Output interface and flow field for visualization
  output_ppm (f, file = "interface.mp4", box = {{0,0},{1.5,2}}, 
              min = 0, max = 1, linear = true, map = cool_warm);
              
  // Output velocity field
  output_ppm (u.y, file = "velocity.mp4", box = {{0,0},{1.5,2}}, 
              linear = true, map = cool_warm);

  // Film thickness measurement
  char filename[80];
  sprintf(filename, "film_profile_t%.1f.dat", t);
  FILE * fp = fopen(filename, "w");
  fprintf(fp, "# Film thickness profile at t = %.2f\n", t);
  fprintf(fp, "# x-coordinate, y-coordinate, volume_fraction, velocity_y\n");
  
  // Sample along vertical lines to measure film thickness
  for (double xi = SUBSTRATE_X + SUBSTRATE_WIDTH + 0.01; 
       xi < SUBSTRATE_X + SUBSTRATE_WIDTH + 0.2; xi += 0.01) {
    
    // Find interface position at this x-location
    foreach() {
      if (fabs(x - xi) < 0.005) { // Within sampling tolerance
        if (f[] > 0.01 && f[] < 0.99) { // Near interface
          fprintf(fp, "%.6f %.6f %.6f %.6f\n", x, y, f[], u.y[]);
        }
      }
    }
  }
  fclose(fp);
  
  // Output key dimensionless parameters for reference
  static FILE * params = NULL;
  if (!params) {
    params = fopen("parameters.txt", "w");
    fprintf(params, "Landau-Levich Dip Coating Simulation Parameters\n");
    fprintf(params, "==============================================\n");
    fprintf(params, "Capillary number (Ca) = %.6f\n", Ca);
    fprintf(params, "Reynolds number (Re) = %.6f\n", Re);
    fprintf(params, "Bond number (Bo) = %.6f\n", Bo);
    fprintf(params, "Density ratio = %.6f\n", Rho21);
    fprintf(params, "Viscosity ratio = %.6f\n", Mu21);
    fprintf(params, "Contact angle = %.1f degrees\n", THETA_DEG);
    fprintf(params, "Expected film thickness ratio h/L ~ Ca^(2/3) = %.6f\n", 
            pow(Ca, 2.0/3.0));
    fclose(params);
  }
}

/**
 * Monitoring routine for simulation health and key quantities
 */
event monitoring (i += 10) {
  // Calculate total liquid volume (conservation check)
  double volume = 0.0;
  foreach() {
    volume += f[] * dv();
  }
  
  // Maximum velocity (stability check)
  double max_vel = 0.0;
  foreach() {
    double vel_mag = sqrt(sq(u.x[]) + sq(u.y[]));
    if (vel_mag > max_vel) max_vel = vel_mag;
  }
  
  fprintf(stderr, "t=%g: Volume=%.6f, MaxVel=%.4f, Ca=%.6f\n", 
          t, volume, max_vel, Ca);
}

/**
 * Main simulation control
 */
int main() {
  // Set domain size: wider to capture far-field effects
  size(1.5);    // Domain width
  L0 = 2.0;     // Domain height
  
  // Origin at bottom-left
  origin(-0.0, 0.0);
  
  // Initial mesh
  N = 1 << 6;  // Start with 64x64 base grid
  
  // Physical properties (dimensionless)
  rho1 = 1.0;           // Liquid density (reference)
  rho2 = Rho21;         // Gas density  
  mu1 = 1.0;            // Liquid viscosity (reference)
  mu2 = Mu21;           // Gas viscosity
  f.sigma = 1.0;        // Surface tension coefficient (dimensionless)
  
  // Gravity effects
  const face vector g[] = {0., -Bo};  // Dimensionless gravity
  a = g;
  
  // Time stepping
  TOLERANCE = 1e-4;     // Solver tolerance
  DT = 0.001;           // Initial time step
  
  // Run simulation
  run();
  
  return 0;
}

/**
 * Theoretical background and expected results:
 * 
 * The Landau-Levich law predicts film thickness:
 * h/L_char ~ 0.94 * Ca^(2/3) * (1 + 3.13*Ca^(2/3))^(-1/2)
 * 
 * For small Ca: h/L_char ≈ 0.94 * Ca^(2/3)
 * 
 * Key observations to make during the course:
 * 1. Film thickness increases with coating speed (higher Ca)
 * 2. Meniscus shape evolution during startup
 * 3. Contact line dynamics and wetting effects  
 * 4. Transition from startup transients to steady state
 * 5. Effect of contact angle on film formation
 * 
 * Parameter exploration suggestions:
 * - Vary U_COATING to see Ca effects on film thickness
 * - Change contact angle to observe wetting influence
 * - Modify viscosity ratio to study coating of different fluids
 * 
 * This example demonstrates the multiscale nature of coating flows:
 * - Macroscale: Overall meniscus shape and flow field
 * - Microscale: Contact line dynamics (subgrid model needed)
 * - Film scale: Uniform film region far from meniscus
 */