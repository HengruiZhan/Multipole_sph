import numpy as np
from scipy.special import sph_harm


class Multipole():
    # def __init__(self, grid, n_moments, dr, center=(0.0, 0.0)):
    def __init__(self, grid, n_moments, dr, center=(0.0, 0.0)):

        self.g = grid
        self.n_moments = n_moments
        self.dr_mp = dr
        self.center = center

        # compute the bins
        # this computation method is not correct
        r_max = max(abs(self.g.rlim[0] - center[0]), abs(self.g.rlim[1] -
                                                         center[0]))
        z_max = max(abs(self.g.zlim[0] - center[1]), abs(self.g.zlim[1] -
                                                         center[1]))

        dmax = np.sqrt(r_max**2 + z_max**2)

        self.n_bins = int(dmax/dr)

        # bin boundaries
        self.r_bin = np.linspace(0.0, dmax, self.n_bins)

        # storage for the inner and outer multipole moment functions
        # we'll index the list by multipole moment l

        self.m_r = []
        self.m_i = []
        for l in range(self.n_moments):
            self.m_r.append(np.zeros((self.n_bins), dtype=np.complex128))
            self.m_i.append(np.zeros((self.n_bins), dtype=np.complex128))

    def compute_expansion(self, rho, l):
        # rho is density that lives on a grid self.g

        # loop over cells
        for i in range(self.g.nr):
            for j in range(self.g.nz):

                # for each cell, i,j, compute r and theta (polar angle from z)
                # and determine which shell we are in
                radius = np.sqrt((self.g.r[i] - self.center[0])**2 +
                                 (self.g.z[j] - self.center[1])**2)

                # tan(theta) = r/z
                theta = np.arctan2(self.g.r[i], self.g.z[j])

                # loop over the multipole moments, l (m = 0 here)
                m_zone = rho[i, j] * self.g.vol[i, j]

                # compute Y_l^m (note: we use theta as the polar
                # angle, scipy is opposite)
                Y_lm = sph_harm(0, l, 0.0, theta)

                R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
                I_lm = np.sqrt(4*np.pi/(2*l + 1)) * Y_lm / radius**(l+1)

                # add to the all of the appropriate inner or outer
                # moment functions
                imask = radius <= self.r_bin
                omask = radius > self.r_bin
                self.m_r[l][imask] += R_lm * m_zone
                self.m_i[l][omask] += I_lm * m_zone

    def sample_mtilde(self, l, r):
        # this returns the result of Eq. 19

        # we need to find which be we are in
        mu_m = np.argwhere(self.r_bin <= r)[-1][0]
        mu_p = np.argwhere(self.r_bin > r)[0][0]

        assert mu_p == mu_m + 1

        mtilde_r = (r - self.r_bin[mu_m])/(self.r_bin[mu_p] - self.r_bin[mu_m]
                                           ) * self.m_r[l][mu_p] + \
            (r - self.r_bin[mu_p])/(self.r_bin[mu_m] -
                                    self.r_bin[mu_p]) * self.m_r[l][mu_m]

        mtilde_i = (r - self.r_bin[mu_m])/(self.r_bin[mu_p] - self.r_bin[mu_m]
                                           ) * self.m_i[l][mu_p] + \
            (r - self.r_bin[mu_p])/(self.r_bin[mu_m] -
                                    self.r_bin[mu_p]) * self.m_i[l][mu_m]

        return mtilde_r, mtilde_i

    def phi(self, r, z, dr, dz=0):
        # return Phi(r), using Eq. 20

        assert dr <= self.g.dr/2 and dr >= 0
        assert dz <= self.g.dz/2 and dz >= 0

        r -= dr

        radius = np.sqrt((r - self.center[0])**2 +
                         (z - self.center[1])**2)

        # tan(theta) = r/z
        theta = np.arctan2(r, z)

        phi_zone = 0.0
        for l in range(self.n_moments):
            mtilde_r, mtilde_i = self.sample_mtilde(l, radius)

            Y_lm = sph_harm(0, l, 0.0, theta)
            R_lm = np.sqrt(4*np.pi/(2*l + 1)) * radius**l * Y_lm
            I_lm = np.nan_to_num(np.sqrt(4*np.pi/(2*l + 1)) *
                                 Y_lm / radius**(l+1))

            phi_zone += sc.G * (mtilde_r * np.conj(I_lm) +
                                np.conj(mtilde_i) * R_lm)

        return -np.real(phi_zone)
