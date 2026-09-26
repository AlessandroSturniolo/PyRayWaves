# PyRayWaves
Simulating optical systems with Python for Physics demonstrations

The project includes two main sub-folders:
 * ```rayoptics```: applies to geometric optics, ray tracing is approached by ray transfer matrix analysis;
 * ```waveoptics```: applies to wave optics (work in progress).

## Ray optics
Sub-module specialised in ray tracing in geometric optical system. The current state-of-the-art implements that by ray transfer matrix analysis, under the paraxial approximation. Each beam is defined by:
 * Initial position (```p```);
 * Initial momentum (```k```);
 * Wavelength (```wavelength```).
From where it can be framed as a vector with components:
 * x-coordinate (```x```);
 * y-coordinate (```y```);
 * Angle between x and beam axis (```theta_x```);
 * Angle between y and beam axis (```theta_y```).

Simple optical components are implemented as matrix transformations, as of now:
 * Thin lenses;
 * Refractive interfaces (planar/spherical);
 * Plane mirror;
 * Thick lens (approximated as a thin lens sandwiched between two spherical refractive interfaces).

## Wave optics
(Under development...)
