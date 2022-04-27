import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import gnlse
from gnlse.common import c

matplotlib.rc('font', size=15)
matplotlib.rc('axes', titlesize=18)
plt.rc('legend', fontsize=15)


if __name__ == '__main__':
    setup = gnlse.GNLSESetup()

    # Physical parameters
    wavelength = 1560  # nm
    total_fiber_length = 21  # m -> total fiber length

    # Input pulse parameters for both modes
    FWHM = 0.023  # ps
    average_power = 0.041  # W
    repetition_rate = 80e6  # Hz
    input_energy = average_power / repetition_rate * 1e9  # nJ

    # Visualization
    ###########################################################################
    theta_list = [1, 89]
    count = len(theta_list)

    plt.figure(figsize=(14, 8), facecolor='w', edgecolor='k')

    for i, theta in enumerate(theta_list):
        # load result
        solution = gnlse.Solution()
        solution.from_file(os.path.join('data',
                                        f'191119_polarisation'
                                        f'_{total_fiber_length}m'
                                        f'_lambda_{wavelength}nm'
                                        f'_power_{int(average_power * 1000)}mW'
                                        f'_angle_{int(theta)}deg_exp.mat'))

        # prepare initial vectors
        Z = solution.Z
        Nz = solution.At.shape[0] // 2
        Atx = solution.At[:Nz, :]
        Aty = solution.At[Nz:, :]
        AWx = solution.AW[:Nz, :]
        AWy = solution.AW[Nz:, :]

        # visualization
        WL_range = [1500, 2200]
        decybl = [-40, 30]
        plt.subplot(3, count, i + 1)
        if i == 0:
            plt.title("(a) Slow axis excitation")
        else:
            plt.title("(b) Fast axis excitation")
        WL = 2 * np.pi * c / solution.W  # wavelength grid
        plt.plot(WL, 10 * np.log10(np.abs(AWx[99, :])**2,
                 where=(np.abs(AWx[99, :])**2 > 0)), '--r',
                 label="2m")
        plt.plot(WL, 10 * np.log10(np.abs(AWx[393, :])**2,
                 where=(np.abs(AWx[393, :])**2 > 0)), '-.r',
                 label="8m")
        plt.plot(WL, 10 * np.log10(np.abs(AWx[981, :])**2,
                 where=(np.abs(AWx[981, :])**2 > 0)), 'r',
                 label="20m")
        plt.plot(WL, 10 * np.log10(np.abs(AWy[99, :])**2,
                 where=(np.abs(AWy[99, :])**2 > 0)), '--g')
        plt.plot(WL, 10 * np.log10(np.abs(AWy[393, :])**2,
                 where=(np.abs(AWy[393, :])**2 > 0)), '-.g')
        plt.plot(WL, 10 * np.log10(np.abs(AWy[981, :])**2,
                 where=(np.abs(AWy[981, :])**2 > 0)), 'g')

        if i == 1:
            plt.legend(fancybox=True, shadow=True)
        else:
            plt.ylabel('Spectral density [dB]')
        plt.xlim(WL_range)
        plt.ylim(decybl)

        plt.subplot(3, count, i + 1 + count)
        solution.AW = AWx
        plt.title("Slow axis")
        gnlse.plot_wavelength_vs_distance(solution, WL_range=WL_range,
                                          norm=1)
        if i == 1:
            plt.ylabel('')
        if ((i + 1 + count) == 3) or ((i + 1 + count) == 4):
            plt.xlabel('')

        plt.subplot(3, count, i + 3 + count)
        solution.AW = AWy
        plt.title("Fast axis")
        gnlse.plot_wavelength_vs_distance(solution, WL_range=WL_range,
                                          norm=1)
        if i == 1:
            plt.ylabel('')

    plt.tight_layout()
    plt.savefig(f'191119_polarisation_{total_fiber_length}m'
                f'_lambda_{wavelength}nm'
                f'_power_{int(average_power * 1000)}mW.png')
    plt.show()
