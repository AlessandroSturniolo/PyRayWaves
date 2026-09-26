import numpy as np

class Ray:

    def __init__(self, k=(0.,0.,1.), p=(0.,0.,0.), wavelength=550.):
        
        if not np.linalg.norm(k) == 1:
            k /= np.linalg.norm(k)

        self.kx = k[0]
        self.ky = k[1]
        self.kz = k[2]
        self.kT = np.sqrt(k[0]**2 + k[1]**2)

        self.r     = np.sqrt(p[0]**2 + p[1]**2)
        self.theta = np.arctan(self.kT / self.kz) if self.kz != 0 else np.pi/2
        self.phi   = np.arctan(self.ky / self.kx) if self.kx != 0 else 0

        self.theta_x = np.arctan(self.kx / self.kz) if self.kz != 0 else np.pi/2
        self.theta_y = np.arctan(self.ky / self.kz) if self.kz != 0 else np.pi/2

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
                 self.kx, self.ky, self.kz
                ]
            ]
        )

        self.wavelength = wavelength
    

    def propagate(self, d):

        M = np.array(
            [[1.0, 0.0,   d, 0.0],
             [0.0, 1.0, 0.0,   d],
             [0.0, 0.0, 1.0, 0.0],
             [0.0, 0.0, 0.0, 1.0]]
        )
        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.z = self.z + (self.kz / np.abs(self.kz))*d
        self.updatePosition()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD
    

    def thinLens(self, f):

        M = np.array(
            [[ 1.0, 0.0, 0.0, 0.0],
             [ 0.0, 1.0, 0.0, 0.0],
             [-1/f, 0.0, 1.0, 0.0],
             [ 0.0,-1/f, 0.0, 1.0]]
        )
        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.updatePosition()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD

    
    def planarInterfaceRefraction(self, n1, n2):

        M = np.array(
            [[ 1.0, 0.0,   0.0,     0.0],
             [ 0.0, 1.0,   0.0,     0.0],
             [ 0.0, 0.0, n1/n2,     0.0],
             [ 0.0, 0.0,   0.0,   n1/n2]]
        )
        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.updatePosition()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )

        return self.rayVectorABCD
    

    def sphericalInterfaceRefraction(self, R, n1, n2):

        a = (n1 - n2)/(n2*R)
        M = np.array(
            [[ 1.0, 0.0,   0.0,   0.0],
             [ 0.0, 1.0,   0.0,   0.0],
             [   a, 0.0, n1/n2,   0.0],
             [ 0.0,   a,   0.0, n1/n2]]
        )

        self.rayVectorABCD = np.dot(M, self.getRayVectorABCD())
        self.updatePosition()
        self.path = np.concatenate(
            (self.path, 
             [
                [
                 self.x, self.y, self.z,
                 self.kx, self.ky, self.kz
                ]
             ] 
            )
        )
        
        return self.rayVectorABCD
    

    def thickLens(self, R1, R2, t, n1, n2):

        self.sphericalInterfaceRefraction(R1,n1,n2)
        self.propagate(t)
        self.sphericalInterfaceRefraction(R2,n2,n1)
        self.updatePosition()

        return self.rayVectorABCD


    def planeReflection(self):

        self.kz = -self.kz
        self.updatePosition()

        return self.rayVectorABCD


    def updatePosition(self):

        self.x, self.y, self.theta_x, self.theta_y = self.rayVectorABCD

        self.kx = self.kz*self.theta_x
        self.ky = self.kz*self.theta_y

        self.kT = np.sqrt(self.kx**2 + self.ky**2)
        
        self.r     = np.sqrt(self.x**2 + self.y**2)
        self.theta = np.arctan(self.kT / self.kz) if self.kz != 0 else np.pi/2
        self.phi   = np.arctan(self.ky / self.kx) if self.kx != 0 else 0


    def getRayVectorABCD(self):
    
        return self.rayVectorABCD


    def getPath(self):

        return self.path
