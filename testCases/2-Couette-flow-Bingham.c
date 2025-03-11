/**
#  Planar Couette flow of Generalized Newtonian Fluid

This code extends the method used in [/sandbox/M1EMN/Exemples/bingham_simple.c](/../sandbox/M1EMN/Exemples/bingham_simple.c)
and generalizes it for any Power Law fluid (using regularization method). Another
difference between the two is that this code calculates the second invariant of
deformation tensor at the face-centers of the cells instead of the cell centers.

## Mathematical Formulations

Unlike the Newtonian fluids, non-Newtonian fluids do not have a linear stress-strain rate relationship.
One way to represent the relationship is using the Generalized Newtonian fluid method:
$$ \tau = \tau_y + 2\mu_0D_{ij}^n $$
The fluid is such that
if $\|\tau\| \le  \tau_y$ then there is no motion $D_{ij}=0\:\forall\:(i,j)$
if the stress is high enough $\|\tau\| >  \tau_y$ then there is motion

*Note:* that $\|\tau\|$ is the modulus defined as the Euclidian norm  $\sqrt{\frac{1}{2}{\tau_{ij} \tau_{ij}}}$.
 It is not $\sqrt{\tau_{11}^2 + \tau_{12}^2}$ as in Balmorth et al. (2006), which is the Frobenius norm.

$D_{ij}$ is the shear strain rate tensor (or the deformation tensor)

$D_{ij}=(u_{i,j}+u_{j,i})/2$: the components in 2D:
$$D_{11}=\frac{\partial u}{\partial x}$$
$$D_{12} =\frac{1}{2}\left(\frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{21} =D_{12} =\frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
$$D_{22}=\frac{\partial v}{\partial y}$$

In the Euclidian norm we have:
$$\|D\|=\sqrt{\frac{D_{ij}D_{ij}}{2}}$$
The second invariant defined by $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
is given by:
$$D_2^2= D_{ij}D_{ij}= \left( \frac{\partial u}{\partial x}\right)^2 + \left(\frac{\partial v}{\partial y}\right)^2
 +  \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)^2$$
and we have obviously $\|D_{ij}\| = D_2/\sqrt{2}$

## Numerical regularization
$$ \tau_{ij} = \tau_y\left(\frac{D_{ij}}{\|D_{ij}\|}\right) + 2\mu_0\|D_{ij}\|^{n-1}D_{ij}^n $$
Factorising with $2D_{ij}$ to obtain a equivalent viscosity
$$\tau_{ij} = 2\left(\mu_0 \|D_{ij}\|^{n-1} + \frac{\tau_y}{2 \|D_{ij}\|}\right)D_{ij}$$
$$\tau_{ij} = 2 \mu_{eq}D_{ij}$$
$$\mu_{eq} = \mu_0\|D_{ij}\|^{n-1} + \frac{\tau_y}{2\|D_{ij}\|}$$
$\mu$ is the min of $\mu_{eq}$ and a large $\mu_{max}$ so that the viscosity does not blow up.
$$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$
*Note:* We present here the formulation in Balmforth, he uses $\dot{\gamma}$ which is by his definition $\sqrt{\frac{1}{2}\dot{\gamma_{ij}}\dot{\gamma_{ij}}}$
and as $\dot{\gamma_{ij}}=2 D_{ij}$ then $\dot{\gamma}$ is $\sqrt{2}D_2$, that is why we have a $\sqrt{2}$ in the equations.


##  Exact solution in the proposed case


We look at an unidirectional flow, a pure shear flow  $u(y)$, $v=0$, so
$D_{11}=D_{22}=0$ and $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$.
$$ \tau_{12} = 2\mu D_{12}^n  + \tau_y =
  2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n + \tau_y $$

Equilibrium between pressure gradient and viscosity (writting $\tau$ for a shorthand of $\tau_{12}$)
$$0=-\frac{\partial p}{\partial x} + \frac{\partial \tau}{\partial y}$$
as there is no stress at the free surface $y=h$, the stress is
$$ \tau = \left(-\frac{\partial p}{\partial x}\right)(h-y)$$
the stress $\tau$ increases from the free surface, as long as $\tau<\tau_y$,
we are under the threshold,
so shear is zero: $\frac{\partial u}{\partial y} =0$,  hence velocity is constant, say it is $U$.
Let us define
 $Y=h-\tau_y/(-\frac{\partial p}{\partial x})$, where  $\tau=\tau_y$.


So :
$$ \left\{\tau<\tau_y, \frac{\partial u}{\partial y} = 0,\:\&\:u=U\:\forall\:Y<y<h\right\} $$
Then going down:
$0<y<Y$ we have $\tau = 2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n + \tau_y$.
This gives:
$$\tau_y + 2^{1-n}\mu\left(\frac{\partial u}{\partial y}\right)^n =  \left(-\frac{\partial p}{\partial x}\right)(Y-y)$$
After some straight forward manipulations:
$$\left(\frac{\partial u}{\partial y}\right) =  \left(\frac{1}{2^{1-n}\mu}\right)^{1/n}\left(-\frac{\partial p}{\partial x}\right)^{1/n}(Y-y)^{1/n}
 = A_{p\mu}(Y-y)^{1/n}$$
and this allows to solve for the velocity profile
$$u = \frac{nA_{p\mu}}{n+1}\left(Y^{\frac{n+1}{n}} - \left(Y-y\right)^{\frac{n+1}{n}}\right)$$
which is indeed zero in $y=0$, and for  $y=Y$, we have the plug flow $u=U$ of value:
$$U= \frac{n}{n+1}A_{p\mu}Y^{\frac{n+1}{n}}$$
$$A_{p\mu} = \left(\frac{1}{2^{1-n}\mu}\right)^{1/n}\left(-\frac{\partial p}{\partial x}\right)^{1/n}$$
For the present case $-\frac{\partial p}{\partial x} = 1, \mu = 1, h = 1$, which gives:
$A_{p\mu} = \frac{1}{2^{(1-n)/n}}, Y = 1 - \tau_y$

# Code
*/
#include "navier-stokes/centered.h"

// Global parameters
char file_name[80];
double tau_y = 0.0;        // Yield stress
double mu_0 = 1.0;         // Base viscosity
double mu_max = 1000.;     // Maximum viscosity for regularization
double n = 1.0;            // Power law exponent
int max_iter = 1e4;        // Maximum iterations
#define DT_MAX (1e-3)      // Maximum timestep

int main(int argc, char const *argv[])
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s filename [mu_0] [tau_y] [n]\n", argv[0]);
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
 Values of yield stress, viscosity, and coefficient.
 - Newtonian: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 1
 - Power law: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 0.5
 - Herschel-Bulkley: $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 0.5
 - Bingham: $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 1
*/
  // Parse command line arguments if provided
  if (argc >= 3) {
    mu_0 = atof(argv[2]);
  }
  if (argc >= 4) {
    tau_y = atof(argv[3]);
  }
  if (argc >= 5) {
    n = atof(argv[4]);
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
// muv will be used as the face vector for viscosity
face vector muv[];

/**
## Initialization event
*/
event init(t = 0) {
  // Set Non-Newtonian viscosity
  mu = muv;
  
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
## Calculating viscosity for generalized Newtonian fluid
*/
event properties(i++) {
  /**
  Implementation of generalized Newtonian viscosity:
  
  The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$ (Frobenius norm)
  $$D_2^2 = D_{ij}D_{ij} = D_{11}^2 + 2D_{12}^2 + D_{22}^2$$
  
  The equivalent viscosity is:
  $$\mu_{eq} = \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{n-1} + 
               \frac{\tau_y}{\sqrt{2} D_2}$$
  
  Finally: $\mu = \min(\mu_{eq}, \mu_{max})$
  */
  foreach_face() {
    // Calculate deformation tensor components at face centers
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1] - u.y[0,-1]) + (u.y[-1,1] - u.y[-1,-1])) / 4.0;
    double D12 = 0.5 * (((u.x[0,1] - u.x[0,-1]) + 
              (u.x[-1,1] - u.x[-1,-1])) / 4.0 + (u.y[] - u.y[-1,0]));
    
    // Calculate second invariant
    double D2 = sqrt(sq(D11) + sq(D22) + 2.0 * sq(D12)) / Delta;
    
    // Calculate effective viscosity
    double mu_temp;
    if (D2 > 0.0) {
      double temp = tau_y / (sqrt(2.0) * D2) + 
                   mu_0 * exp((n - 1.0) * log(D2 / sqrt(2.0)));
      mu_temp = min(temp, mu_max);
    } else {
      if (tau_y > 0.0 || n < 1.0) {
        mu_temp = mu_max;
      } else {
        mu_temp = (n == 1.0 ? mu_0 : 0.0);
      }
    }
    
    // Apply viscosity at face
    muv.x[] = fm.x[] * mu_temp;
  }
}