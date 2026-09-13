import numpy as np
import cmath

class Wave:

    def __init__(self, E0, phase, k=(0.,0.,1.), wavelength=550.):
        self.E0    = np.array(E0)
        self.phase = np.array(phase)
        self.J     = E0*np.exp(1.0j*self.phase)  # Jones' vector (complex)
        self.I     = np.sum(self.E0**2)

        if not np.linalg.norm(k) == 1:
            k /= np.linalg.norm(k)

        self.k = k*2*np.pi/wavelength
        self.c = 3.0e+17 # Speed of light [nm/s]
        self.omega = self.c*2*np.pi/wavelength

    def propagate(self, x, t):

        # x is expected to be an array-like object
        x = np.array(x)
        return self.J*np.exp(1.0j*(np.dot(self.k, x) + self.omega*t))
