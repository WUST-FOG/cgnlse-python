"""
Example of modulation instability in highly birefringent
fibers with circularly polarized eigenmodes [1].

[1] Zołnacz, K., Tarnowski, K. L., Napiórkowski, M.,
    Poturaj, K., Mergo, P., Urbanczyk, W.
    Vector modulation instability in highly birefringent
    fibers with circularly polarized eigenmodes.
    IEEE Photonics Journal, 13(1):7100616, 2021.
"""

import numpy as np
import matplotlib.pyplot as plt

import gnlse

from cnlse import CNLSE, DubbleDispersionFiberFromTaylor

if __name__ == '__main__':
    setup = gnlse.GNLSESetup()

    # Numerical parameters
    setup.resolution = 2**11
    setup.time_window = 20  # ps
    setup.z_saves = 200

    # Physical parameters
    setup.wavelength = 1064  # nm
    w0 = 2 * np.pi * gnlse.common.c / setup.wavelength  # 1/ps = THz
    setup.fiber_length = 2  # m
    n2 = 2.6e-20  # m^2/W
    Aeff0 = 54e-12  # 1/m^2
    gamma = n2 * 2 * np.pi / setup.wavelength * 1e9 / Aeff0  # 1/W/m
    setup.nonlinearity = (gamma, gamma)

    # The dispersion model is built from a Taylor expansion
    # with coefficients given below.
    loss = 0
    D = -40  # ps/km/nm
    beta2 = -setup.wavelength**2 * D / 2 / np.pi / gnlse.common.c * 1e-3
    betas = np.array([[beta2], [beta2]])
    G = 5e-4
    Dbeta1 = G / (gnlse.common.c / 1e9)  # ps/m
    dn = 1e-4
    setup.dispersion_model = DubbleDispersionFiberFromTaylor(
        loss, betas, dn, Dbeta1, w0)

    # Input pulse parameters
    peak_power = 2000  # W
    peak_noise = 1e-7  # W
    # Time domain grid
    t = np.linspace(-setup.time_window / 2,
                    setup.time_window / 2,
                    setup.resolution)
    At = gnlse.CWEnvelope(
        peak_power / 2, peak_noise / 2).A(np.concatenate((t, t)))
    setup.pulse_model = At

    # Nonlinear Simulation
    ###########################################################################
    solver = CNLSE(setup)
    solver.report = True
    solution = solver.run()

    # Visualization
    ###########################################################################
    # prepare initial vectors
    Z = solution.Z
    Atx = solution.At[:setup.z_saves, :]
    Aty = solution.At[setup.z_saves:, :]
    AWx = solution.AW[:setup.z_saves, :]
    AWy = solution.AW[setup.z_saves:, :]

    lambdas = 2 * np.pi * gnlse.common.c / solution.W
    V = solution.W - w0
    gx = np.log(abs(AWx[-1, :]
                    )**2 / abs(AWx[0, :])**2) / setup.fiber_length
    gy = np.log(abs(AWy[-1, :]
                    )**2 / abs(AWy[0, :])**2) / setup.fiber_length

    plt.figure(figsize=(10, 10))

    plt.subplot(3, 1, 1)
    plt.title("Results for modulation instability")
    plt.plot(solution.t, abs(Atx[0, :])**2)
    plt.plot(solution.t, abs(Aty[0, :])**2)
    plt.xlabel('t [ps]')
    plt.ylabel(r'$|A|^2$ [W]')

    plt.subplot(3, 1, 2)
    plt.semilogy(lambdas, abs(AWx[-1, :])**2)
    plt.semilogy(lambdas, abs(AWy[-1, :])**2)
    plt.xlim([950, 1200])
    plt.xlabel(r'$\lambda$ [nm]')
    plt.ylabel('Intensity [a.u.]')

    plt.subplot(3, 1, 3)
    plt.plot(V, gx)
    plt.plot(V, gy)
    plt.xlabel(r'$\Omega$ [THz]')
    plt.ylabel('Gain [1/m]')

    plt.show()
