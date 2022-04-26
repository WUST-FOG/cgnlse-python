import numpy as np
import scipy.integrate
from scipy.fftpack import fft, ifft, fftshift
import tqdm

from gnlse import GNLSESetup, Solution
from gnlse.common import c


class CNLSE:
    """
    Model propagation of an optical pulse in a fiber by integrating
    the two coupled generalized non-linear Schrödinger equations.

    Attributes
    ----------
    setup : GNLSESetup
        Model inputs in the form of a ``GNLSESetup`` object.
    """
    def __init__(self, setup):
        if not isinstance(setup, GNLSESetup):
            raise TypeError("setup is not an instance of GNLSESetup")

        if setup.resolution is None:
            raise ValueError("'resolution' not set")
        if setup.time_window is None:
            raise ValueError("'time_window' not set")
        if setup.wavelength is None:
            raise ValueError("'wavelength' not set")
        if setup.fiber_length is None:
            raise ValueError("'fiber_length' not set")
        if setup.pulse_model is None:
            raise ValueError("'pulse_model' not set")

        # simulation parameters
        self.fiber_length = setup.fiber_length
        self.z_saves = setup.z_saves
        self.rtol = setup.rtol
        self.atol = setup.atol
        self.method = setup.method
        self.N = setup.resolution
        self.report = False

        # Time domain grid
        self.t = np.linspace(-setup.time_window / 2,
                             setup.time_window / 2,
                             self.N)

        # Relative angular frequency grid
        self.V = 2 * np.pi * np.arange(-self.N / 2,
                                       self.N / 2
                                       ) / (self.N * (self.t[1] - self.t[0]))
        # Central angular frequency [10^12 rad]
        w_0 = (2.0 * np.pi * c) / setup.wavelength
        self.Omega = self.V + w_0

        # Absolute angular frequency grid
        self.W = fftshift(self.Omega)

        # Raman scattering
        self.fr, h1T, h2T, h3T = setup.raman_model(self.t)
        self.h1W = self.N * ifft(
            fftshift(h1T))
        self.h2W = self.N * ifft(
            fftshift(h2T))
        self.h3W = self.N * ifft(
            fftshift(h3T))

        # Input pulse
        if hasattr(setup.pulse_model, 'A'):
            self.A = setup.pulse_model.A(self.t)
        else:
            self.A = setup.pulse_model

        # Dispersion operator
        self.D, self.deltaB = setup.dispersion_model.D(self.V)

        # Nonlinearity
        gamma, scale = setup.nonlinearity[0].gamma(self.V)
        self.gammax = fftshift(gamma / w_0)
        self.scale = fftshift(scale)
        gamma, _ = setup.nonlinearity[1].gamma(self.V)
        self.gammay = fftshift(gamma / w_0)

    def run(self):
        """
        Solve two mode CNLSE equation.

        Returns
        -------
        setup : Solution
            Simulation results in the form of a ``Solution`` object.
        """
        if self.A.size != self.N * 2:
            raise ValueError("'pulse_model' has not enougth values")

        if self.D.size != self.N * 2:
            raise ValueError("'dispersion' has not enougth values")

        self.Dx = fftshift(self.D[:self.N])
        self.Dy = fftshift(self.D[self.N:])
        dt = self.t[1] - self.t[0]
        Z = np.linspace(0, self.fiber_length, self.z_saves)

        if self.report:
            progress_bar = tqdm.tqdm(total=self.fiber_length, unit='m')

        def rhs(z, AW):
            """
            The right hand side of the differential equation to integrate.
            """
            if self.report:
                progress_bar.n = round(z, 3)
                progress_bar.update(0)

            CxW = AW[:self.N]
            CyW = AW[self.N:]
            Atx = fft(CxW * np.exp(self.Dx * z))
            Aty = fft(CyW * np.exp(self.Dy * z))
            ITx = np.abs(Atx)**2
            ITy = np.abs(Aty)**2

            exp_m2i_deltaB_z = np.exp(-2j * self.deltaB * z)
            exp_m1i_deltaB_z = np.exp(-1j * self.deltaB * z)
            exp_p1i_deltaB_z = np.exp(+1j * self.deltaB * z)
            exp_p2i_deltaB_z = np.exp(+2j * self.deltaB * z)

            Mx = ifft((1 - self.fr) * ((ITx + 2 / 3 * ITy)
                    * Atx + 1 / 3 * Aty**2 * np.conj(Atx) * exp_m2i_deltaB_z)
                    + self.fr * dt * (Atx * (fft(self.h1W * ifft(ITx)) +
                    fft(self.h2W * ifft(ITy))
                    ) + Aty * (fft(self.h3W * ifft(
                    Atx * np.conj(Aty) * exp_p1i_deltaB_z
                    + Aty * np.conj(Atx) * exp_m1i_deltaB_z))
                    ) * exp_m1i_deltaB_z))

            My = ifft((1 - self.fr)
                    * (ITy * Aty + 2 / 3 * ITx * Aty
                    + 1 / 3 * Atx**2 * np.conj(Aty) * exp_p2i_deltaB_z)
                    + self.fr * dt * (
                        Aty * (
                    fft(self.h1W * ifft(ITy)) + fft(self.h2W * ifft(ITx))
                ) + Atx * (
                    fft(
                        self.h3W * ifft(
                            Atx * np.conj(Aty) * exp_p1i_deltaB_z
                            + Aty * np.conj(Atx) * exp_m1i_deltaB_z))
                ) * exp_p1i_deltaB_z))

            rx = 1j * self.gammax * self.W * Mx * np.exp(-self.Dx * z)
            ry = 1j * self.gammay * self.W * My * np.exp(-self.Dy * z)

            return np.concatenate((rx, ry))

        Ax = self.A[:self.N]
        Ay = self.A[self.N:]
        solution = scipy.integrate.solve_ivp(
            rhs,
            t_span=(0, self.fiber_length),
            y0=np.concatenate(
                (ifft(Ax) * self.scale,
                 ifft(Ay) * self.scale)),
            t_eval=Z,
            rtol=self.rtol,
            atol=self.atol,
            method=self.method,
            dense_output=True)

        if self.report:
            progress_bar.close()

        AWx = solution.y.T[:, :self.N]
        AWy = solution.y.T[:, self.N:]

        # Transform the results into the time domain
        Atx = np.zeros(AWx.shape, dtype=AWx.dtype)
        Aty = np.zeros(AWy.shape, dtype=AWy.dtype)
        for i in range(len(AWx[:, 0])):
            AWx[i, :] *= np.exp(np.transpose(
                self.Dx) * Z[i]) / self.scale
            Atx[i, :] = fft(AWx[i, :])
            AWx[i, :] = fftshift(AWx[i, :]) * self.N * dt
            AWy[i, :] *= np.exp(np.transpose(
                self.Dy) * Z[i]) / self.scale
            Aty[i, :] = fft(AWy[i, :])
            AWy[i, :] = fftshift(AWy[i, :]) * self.N * dt

        At = np.concatenate((Atx, Aty))
        AW = np.concatenate((AWx, AWy))

        return Solution(self.t, self.Omega, Z, At, AW)
