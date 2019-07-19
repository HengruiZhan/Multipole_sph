import analytical_potential as ap
import numpy as np
import matplotlib.pyplot as plt
import grid
import multipole
from scipy import interpolate
"""
nr = 128
nz = 256
"""
nr = 256
nz = 256

g = grid.Grid(nr, nz)
rlim = (0, 0.1)
zlim = (-1.0, 1.0)
dens = g.scratch_array()
"""
# density of a perfect sphere
sph_center = (0.0, 0.0)
radius = np.sqrt((g.r2d - sph_center[0])**2 + (g.z2d - sph_center[1])**2)
dens[radius <= 0.3] = 1.0
"""

# density of a MacLaurin spheroid
sph_center = (0.0, 0.0)
a_1 = 0.23
a_3 = 0.10
mask = g.r2d**2/a_1**2 + g.z2d**2/a_3**2 <= 1
dens[mask] = 1.0
density = 1.0

phi_mac = g.scratch_array()
for i in range(g.nr):
    for j in range(g.nz):
        mac_phi = ap.Ana_Mac_pot(a_1, a_3, g.r[i], g.z[j], density)
        phi_mac[i, j] = mac_phi.potential
plt.imshow(np.log10(np.abs(np.transpose(phi_mac))), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("phi.png")

"""
plt.imshow(np.transpose(dens), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])

ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("dens.png")
e = np.sqrt(1 - (a_3/a_1)**2)
"""
"""
center = (0.0, 0.0)
n_moments = 40
m = multipole.Multipole(g, n_moments, 1.3*g.dr, center=center)
for l in range(n_moments):
    m.compute_expansion(dens, l)

phi = g.scratch_array()

for i in range(g.nr):
    for j in range(g.nz):
        phi[i, j] = m.phi(g.r[i], g.z[j])

plt.imshow(np.log10(np.abs(np.transpose(phi))), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("phi.png")
"""
# PRACTICE interpolate
# x = (0, 10, 20, 30, 40, 50, 60, 100, 150, 200, 300, 400)
"""
x = (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200)
y = (10, 9.5, 7.2, 4.0, 3.5, 3.3, 3.2, 3.1, 2.9, 2.7, 2.5, 2.3)
# scipy.interpolate.LSQUnivariateSpline is good metod
f = interpolate.interp1d(x, y, kind='cubic')
xnew = np.arange(0, 250, 50)
ynew = f(xnew)
plt.plot(x, y, 'o', xnew, ynew, '-')
plt.show()
"""
