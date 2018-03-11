import scipy as sp
from scipy import special
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.colors import LogNorm
import pylab

# Quantum Numbers
qn = [14]
ql = [8]
qm = [0]

# bohr radius
a = 0.529E-10

# points
density = 100
r_line = np.linspace(0, 6E-9, density, endpoint=True)
theta_line = np.linspace(0, 2*np.pi, density)
r, theta = np.meshgrid(r_line, theta_line)
phi = 0

#fig = plt.figure()

for n in qn:
    for l in ql:
        for m in qm:
            if n<=l:
                continue
            if m>l:
                continue
#            ax = fig.add_subplot(1, 1, 1, projection='3d')
            R = np.sqrt((2/(n*a))*(np.math.factorial(n-l-1)/(2*n*np.math.factorial(n+1))**3))*np.exp(-r/(n*a)) * ( 2*r / (n*a))**l * special.assoc_laguerre( 2*r / (n*a),n-l-1, 2*l+1 )
            ep = (-1)**m
            Y = ep*np.sqrt((2*l+1)/(4*np.pi) * special.factorial(l-np.abs(m))/special.factorial(1+np.abs(m))) * np.exp(1j*m*phi) * special.lpmv(m,l,np.cos(theta))

            Z = np.abs(R*Y)**2
            X, Y = r*np.cos(theta), r*np.sin(theta)

            pylab.pcolor(Y, X, Z, cmap='jet')
            pylab.colorbar()



pylab.show()
