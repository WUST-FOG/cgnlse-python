"""
Calculates Raman responses
"""

import numpy as np


def raman_polarisation(T):
    """Raman scattering function for silica optical fibers, based on
    Prannay Balla and Govind P. Agrawal
    Phys. Rev. A 98, 023822 (2018)
    (J. Opt. Soc. Am. B 6, 1159–1166 (1989),
    J. Opt. Soc. Am. B 9, 1061–1082 (1992))

    Parameters
     ----------
    T : float
        Time vector.

    Returns
    -------
    fr : float
       Share of Raman response.
    h1T : ndarray
       Vector representing h1 temporal Raman response.
    h2T : ndarray
       Vector representing h2 temporal Raman response.
    h3T : ndarray
       Vector representing h3 temporal Raman response.

    """

    # Raman response [arbitrary units]
    fr = 0.245
    # Adjustable parameters used to fit the actual Raman gain spectrum [ps]
    tau1 = 0.0122
    tau2 = 0.032
    taub = 0.096
    # Fractional contribution of the anisotropic reponse to the total Raman
    # response
    fb = 0.21
    fc = 0.04
    # Fractional contribution of the isotropic reponse to the total Raman
    # response
    fa = 1 - fb - fc
    # Anisotropic Raman response
    ha = (tau1**2 + tau2**2) / tau1 / (tau2**2) * np.exp(-T / tau2) * np.sin(
        T / tau1)
    # Izotropic Raman respons
    hb = (2 * taub - T) / (taub**2) * np.exp(-T / taub)

    h1T = (fa + fc) * ha + fb * hb
    h2T = fa * ha
    h3T = (fc * ha + fb * hb) / 2

    h1T[T < 0] = 0
    h2T[T < 0] = 0
    h3T[T < 0] = 0

    return fr, h1T, h2T, h3T
