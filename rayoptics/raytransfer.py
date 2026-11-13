import numpy as np

# ============= RAY TRANSFER MATRIX ANALYSIS ============= #
# The following class implements a simple ray tracing using a ray transfer matrix analysis approach.
# A basic ray object is defined by its starting position and direction, and is handled as a vector in a 4D space (x, y, theta_x, theta_y), 
# where x and y are the transverse coordinates, while theta_x and theta_y are the angles 
# of the ray with respect to the optical axis in the x and y directions, respectively.
#
# Define a beam by:
# - Starting position: p = (x0, y0, z0);
# - Direction versor: k = (kx, ky, kz);
# - Wavelength (in nm, default 550 nm). 
#
# Refractive indices are assumed to be isotropic.
#
# Author: Alessandro Sturniolo

class Ray:

    def __init__(self, k=(0.,0.,1.), p=(0.,0.,0.), wavelength=550.):

        # Starting definition of the ray object

        self.kx = k[0]
        self.ky = k[1]
        self.kz = k[2]

        if self.kz == 0:
            raise ValueError("The z-component of the direction vector must be non-zero.")

        self.theta_x = self.kx / self.kz
        self.theta_y = self.ky / self.kz

        self.x = p[0]
        self.y = p[1]
        self.z = p[2]

        self.rayVectorABCD = np.array(
            [
             self.x, 
             self.y, 
             self.theta_x,
             self.theta_y
            ]
        )

        self.path = np.array(
            [
                [ 
                 self.x,  self.y,  self.z,
                 self.theta_x, self.theta_y,
                 self.kx, self.ky, self.kz
                ]
            ]
        )

        self.wavelength = wavelength
        self.flip_z     = 1.   # Will be inverted in case of reflection
    

    def propagate(self, d):

        # Unperturbed propagation of the ray over a distance d
        M = np.array(
            [[1.0, 0.0,   d, 0.0],
             [0.0, 1.0, 0.0,   d],
             [0.0, 0.0, 1.0, 0.0],
             [0.0, 0.0, 0.0, 1.0]]
        )
        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.z = self.z + self.flip_z*d
        self.updatePositionDirection()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.theta_x, self.theta_y,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD
    

    def thinLens(self, f):

        # Refraction of the ray through a thin lens with focal length f
        M = np.array(
            [[ 1.0, 0.0, 0.0, 0.0],
             [ 0.0, 1.0, 0.0, 0.0],
             [-1/f, 0.0, 1.0, 0.0],
             [ 0.0,-1/f, 0.0, 1.0]]
        )
        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.updatePositionDirection()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.theta_x, self.theta_y,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD

    
    def planarInterfaceRefraction(self, n1, n2):

        # Refraction of the ray through a planar interface (xy-plane) between two different media (refractive indices: n1, n2)
        M = np.array(
            [[ 1.0, 0.0,   0.0,     0.0],
             [ 0.0, 1.0,   0.0,     0.0],
             [ 0.0, 0.0, n1/n2,     0.0],
             [ 0.0, 0.0,   0.0,   n1/n2]]
        )
        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.updatePositionDirection()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.theta_x, self.theta_y,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD
    

    def sphericalInterfaceRefraction(self, R, n1, n2):

        # Refraction of the ray through a spherical interface between two different media (refractive indices: n1, n2), with curvature radius R
        a = (n1 - n2)/(n2*R)
        M = np.array(
            [[ 1.0, 0.0,   0.0,   0.0],
             [ 0.0, 1.0,   0.0,   0.0],
             [   a, 0.0, n1/n2,   0.0],
             [ 0.0,   a,   0.0, n1/n2]]
        )

        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.updatePositionDirection()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.theta_x, self.theta_y,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )
        
        return self.rayVectorABCD
    

    def thickLens(self, R1, R2, t, n1, n2):

        # Propagation through a lens of finite thickness t, with two spherical interfaces of radii R1 and R2
        # Combining refraction through spherical and thick lenses
        self.sphericalInterfaceRefraction(R1,n1,n2)
        self.propagate(t)
        self.sphericalInterfaceRefraction(R2,n2,n1)
        self.updatePositionDirection()

        return self.rayVectorABCD


    def planarMirrorReflection(self):

        # Reflection on a planar mirror (xy-plane)
        self.flip_z = -self.flip_z
        self.updatePositionDirection()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.theta_x, self.theta_y,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD


    def updatePositionDirection(self):
        
        # Updating ray positions and direction
        self.x, self.y, self.theta_x, self.theta_y = self.rayVectorABCD

        self.kx = self.kz*self.theta_x
        self.ky = self.kz*self.theta_y


    def getRayVectorABCD(self):

        # Retrieving the ray vector in the ABCD formalism
        return self.rayVectorABCD


    def getPath(self):

        # Retrieving the full path of the ray
        return self.path
