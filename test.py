# from gekko import GEKKO
import analytical_potential as ap
import comparitive_difference as comdf
import numpy as np
import matplotlib.pyplot as plt
import grid
import L2_difference as L2
import multipole

nr = 256
nz = 512

g = grid.Grid(nr, nz)
rlim = (0, 0.1)
zlim = (-1.0, 1.0)
dens = g.scratch_array()


# density of a perfect sphere
sph_center = (0.0, 0.0)
radius = np.sqrt((g.r2d - sph_center[0])**2 + (g.z2d - sph_center[1])**2)
a = 0.3
dens[radius <= a] = 1.0
# analytical potential of a perfect sphere
density = 1.0
phi_sph = g.scratch_array()
for i in range(g.nr):
    for j in range(g.nz):
        sph_phi = ap.Ana_Sph_pot(g.r[i], g.z[j], a, density)
        phi_sph[i, j] = sph_phi.potential

"""
# density of a MacLaurin spheroid
sph_center = (0.0, 0.0)
a_1 = 0.23
a_3 = 0.10
mask = g.r2d**2/a_1**2 + g.z2d**2/a_3**2 <= 1
dens[mask] = 1.0
density = 1.0
# analytical potential of a MacLaurin spheroid
phi_mac = g.scratch_array()
for i in range(g.nr):
    for j in range(g.nz):
        mac_phi = ap.Ana_Mac_pot(a_1, a_3, g.r[i], g.z[j], density)
        phi_mac[i, j] = mac_phi.potential
"""
"""
# plot the density
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
plt.imshow(np.log10(np.abs(np.transpose(phi_mac))), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])

plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("anal_phi.png")
"""

# create a multipole object
center = (0.0, 0.0)
n_moments = 1
m = multipole.Multipole(g, n_moments, 1.2*g.dr, center=center)
# calculate the multipole
for l in range(n_moments):
    m.compute_expansion(dens, l)

"""
phi = g.scratch_array()
for i in range(g.nr):
    for j in range(g.nz):
        phi[i, j] = m.phi(g.r[i], g.z[j])
"""
# evaluate eq.20 at the surface
# the surface coodinates satisfy:
# r_surface = r_grid_center - dr/2

phi = g.scratch_array()
for i in range(g.nr):
    for j in range(g.nz):
        phi[i, j] = m.phi(g.r[i], g.z[j], 0)
        # phi[i, j] = m.phi(g.r[i], g.z[j], g.dr/2)
"""
plt.imshow(np.log10(np.abs(np.transpose(phi))), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])
plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("lmax=0.png")

L2normerr = L2.L2_diff(phi_mac, phi, g.scratch_array())
print("lmax =", n_moments-1)
print(L2normerr)
"""
# calculate and plot the comparitive difference
diff = comdf.comp_diff(phi_sph, phi)
plt.imshow(np.abs(np.transpose(diff.difference)), origin="lower",
           interpolation="nearest",
           extent=[g.rlim[0], g.rlim[1],
                   g.zlim[0], g.zlim[1]])
plt.colorbar()
ax = plt.gca()
ax.set_aspect("equal")
plt.savefig("comparitive_difference_lmax=0.png")

"""
# linear interpolate the L2_norm error
xp = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
               15, 16, 17, 18, 19, 20])
yp = np.array([1.2779725071002864e-12, 1.2779725071002864e-12,
               1.341920758799657e-12, 1.341920758799657e-12,
               1.3349025961125413e-12, 1.3349025961125413e-12,
               1.3366412358041248e-12, 1.3366412358041248e-12,
               1.3363470028306473e-12, 1.3363470028306473e-12,
               1.3364726905325435e-12, 1.3364726905325435e-12,
               1.3364595642792355e-12, 1.3364595642792355e-12,
               1.3364739556107749e-12, 1.3364739556107749e-12,
               1.3364769845856654e-12, 1.3364769845856654e-12,
               1.3364822456081226e-12, 1.3364822456081226e-12,
               1.33648718522299e-12])
x = np.arange(0, 22)
y = np.interp(x, xp, yp)
plt.plot(xp, yp, 'o')
plt.plot(x, y, 'r--')
plt.show()
"""
