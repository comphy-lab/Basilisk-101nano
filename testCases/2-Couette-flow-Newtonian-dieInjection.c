#include "navier-stokes/centered.h"
char filename[80];
int imax = 1e4;
#define dtmax (1e-3)

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  init_grid (1<<6);
  L0 = 1.0;
  origin (0.0, 0.0);
  DT = dtmax;
  stokes = true;
  TOLERANCE = 1e-5;

/**
 Values of yeild stress, viscosity, and coefficient.<br/>
 Newtonian: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 1<br/>
 Power law $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 0.5<br/>
 Herschel-Bulkley $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 0.5<br/>
 Bingham $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 1<br/>
*/

  periodic (right);
/**
  slip at the top
*/
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  uf.n[top] = neumann(0);
/**
 no slip at the bottom
*/
  u.n[bottom] = dirichlet(0);
  uf.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
/**
 presure conditions are neumann 0.0
 */
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
event init (t = 0) {
  // preparing viscosity to be used as Non-Newtonian fluid
  mu = muv;
  /**
    presure gradient `mdpdx`
   $$-\frac{dp}{dx} = 1 $$
  */
  const face vector mdpdx[] = {1.0,0.0};
  a = mdpdx;
  /**
   Initialy at rest
  */
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
  }
  foreach(){
    un[] = u.x[];
  }
  dump (file = "start");
}

/**
We look for a stationary solution. */
event logfile (i += 500; i <= imax) {
  double du = change (u.x, un);
  fprintf(ferr, "i = %d: err = %g\n", i, du);
  if (i > 0 && du < 1e-6){
    dump (file = filename);
    return 1; /* stop */
  }
  if (i==imax){
    dump (file = filename);
  }
}


event properties(i++) { // Overloading the properties event
  /**
  ## Implementation of generalized Newtonian viscosity

   $$D_{11} = \frac{\partial u}{\partial x}$$
   $$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{22} = \frac{\partial v}{\partial y}$$
   The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$ (this is the Frobenius norm)
   $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
   the equivalent viscosity is
   $$\mu_{eq}= \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{N-1} + \frac{\tau_y}{\sqrt{2} D_2 }$$
   **Note:** $\|D\| = D_2/\sqrt{2}$

   Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$.

   The fluid flows always, it is not a solid, but a very viscous fluid.
   $$ \mu = \text{min}\left(\mu_{eq}, \mu_{max}\right) $$
  */
  double muTemp = mu_0;
  foreach_face() {
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
    double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    if (D2 > 0.0) {
      double temp = tauy/(sqrt(2.0)*D2) + mu_0*exp((n-1.0)*log(D2/sqrt(2.0)));
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.0 || n < 1.0){
        muTemp = mumax;
      } else {
        muTemp = (n == 1.0 ? mu_0 : 0.0);
      }
    }
    muv.x[] = fm.x[]*(muTemp);
  }
}