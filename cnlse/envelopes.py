"""Amplitude envelopes of different pulses

This module contains functions drawing envelopes of various pulses:
hyperbolic secant, gaussian and lorentzian.

"""

import numpy as np
from scipy.fftpack import fft, ifft

from gnlse.common import c, hbar
from gnlse.envelopes import SechEnvelope, Envelope


class DoubleSechEnvelope(Envelope):
    """Amplitude envelope of hyperbolic secant pulse.

    Attributes
    ----------
    input_energy : float
        Input pulse energy [nJ]
    FWHM : float
        Pulse duration Full-Width Half-Maximum [ps]
    theta : float
        Azimuth angle [rad]
    w0 : float
        Central frequency of the input pulse [THz]
    """

    def __init__(self, input_energy, FWHM, theta, w0):
        self.name = 'Double Hyperbolic secant envelope'
        self.energy = input_energy
        self.FWHM = FWHM
        self.theta = theta
        self.w0 = w0

    def A(self, T):
        """

        Parameters
        ----------
        T : ndarray, (n, )
            Time vector

        Returns
        -------
        ndarray, (n, )
            Amplitude envelope of double hyperbolic secant pulse in time.
        """
        N = len(T)
        dT = T[1] - T[0]
        V = 2 * np.pi * np.array(np.arange(-N/2, N/2)) / (N * dT)
        ATx = np.cos(self.theta) * SechEnvelope(1, self.FWHM).A(T)
        ATy = np.sin(self.theta) * SechEnvelope(1, self.FWHM).A(T)
        energy_x = np.trapz(-T, abs(ATx)**2) / 1e3  # nJ
        energy_y = np.trapz(-T, abs(ATy)**2) / 1e3  # nJ
        energy = energy_x + energy_y
        ATx = ATx * np.sqrt(self.energy / energy)
        ATy = ATy * np.sqrt(self.energy / energy)
        AWx = ifft(ATx) * dT * N
        AWy = ifft(ATy) * dT * N
        photon_energy = hbar * (
            V + (2.0 * np.pi * c) / self.w0)  # [J.s.THz = TJ]
        photon_energy = photon_energy * 1e24  # [pJ]
        noise = np.sqrt(photon_energy * dT * N)  # [sqrt(pJ.ps)]
        AWx_noise = noise * np.exp(2j * np.pi * np.random.rand(1, N))
        AWy_noise = noise * np.exp(2j * np.pi * np.random.rand(1, N))
        AWx = AWx + AWx_noise[0, :]
        AWy = AWy + AWy_noise[0, :]
        ATx = fft(AWx / dT / N)
        ATy = fft(AWy / dT / N)
        return np.concatenate((ATx, ATy))
