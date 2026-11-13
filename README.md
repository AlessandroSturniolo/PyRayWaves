# PyRayWaves
Simulating optical systems with Python for Physics demonstrations

The project includes two main sub-folders:
 * ```rayoptics```: reference for geometric optics:
   * Implementing ray tracing with a simple ray transfer matrix analysis approach;
   * Planning to implement Gaussian Beam objects (future development);
 * ```waveoptics```: reference for wave optics:
   * Planning to implement polarisation effects, e.g. Jones calculus (future development).

## Ray tracing
Sub-module specialised in ray tracing in geometric optical system. The current state-of-the-art implements that by ray transfer matrix analysis, under the paraxial approximation. Each beam is defined by:
 * Initial position (```p```);
 * Initial momentum (```k```);
 * Wavelength (```wavelength```, in nanometers).
From where it can be framed as a vector with components:
 * x-coordinate (```x```);
 * y-coordinate (```y```);
 * Angle between x and beam axis (```theta_x```);
 * Angle between y and beam axis (```theta_y```).

Simple optical components are implemented as matrix transformations, as of now:
 * Thin lenses;
 * Refractive interfaces (planar/spherical);
 * Plane mirrors;
 * Thick lens (approximated as a thin lens sandwiched between two spherical refractive interfaces).

## Wave optics
(Under development...)
**Planned:** implementing EM fields, interference effects, Jones calculus.
