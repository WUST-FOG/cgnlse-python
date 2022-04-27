import os
import numpy as np
import tqdm

import gnlse

from cnlse import (CNLSE, DubbleDispersionFiberFromOPD,
                   DoubleSechEnvelope, raman_polarisation)


if __name__ == '__main__':
    setup = gnlse.GNLSESetup()

    # Numerical parameters
    setup.resolution = 2**13
    setup.time_window = 45  # ps
    setup.z_saves = 50
    setup.rtol = 1e-4
    setup.atol = 1e-6
    setup.method = 'RK45'

    # Physical parameters
    setup.wavelength = 1560  # nm
    w0 = (2.0 * np.pi * gnlse.common.c) / setup.wavelength  # THz
    setup.fiber_length = 1  # m -> fiber length to filter temporal noise
    total_fiber_length = 21  # m -> total fiber length
    setup.raman_model = raman_polarisation

    # Input pulse parameters for both modes
    FWHM = 0.023  # ps
    average_power = 0.041  # W
    repetition_rate = 80e6  # Hz
    input_energy = average_power / repetition_rate * 1e9  # nJ

    # Read input files for calculating the nonlinearity
    # read mat file for slow mode
    mat_path = os.path.join(os.path.dirname(__file__),
                            'data', '191119_Eslow_KG.mat')
    mat = gnlse.read_mat(mat_path)
    # neffs
    neffx = mat['neff_r'][:, 0]
    # wavelengths in nm
    lambdasx = mat['lams'][:, 0] * 1e3
    # efective mode area in m^2
    Aeffx = mat['Aeff'][:, 0] * 1e-12

    # read mat file for fast mode
    mat_path = os.path.join(os.path.dirname(__file__),
                            'data', '191119_Efast_KG.mat')
    mat = gnlse.read_mat(mat_path)
    # neffs
    neffy = mat['neff_r'][:, 0]
    # wavelengths in nm
    lambdasy = mat['lams'][:, 0] * 1e3
    # efective mode area in m^2
    Aeffy = mat['Aeff'][:, 0] * 1e-12
    # nonlinear index of refraction [m^2/W]
    n2 = 2.6e-20
    fun_gammax = gnlse.NonlinearityFromEffectiveArea(neffx,
                                                     Aeffx,
                                                     lambdasx,
                                                     setup.wavelength,
                                                     n2,
                                                     neff_max=10)
    fun_gammay = gnlse.NonlinearityFromEffectiveArea(neffy,
                                                     Aeffx,
                                                     lambdasy,
                                                     setup.wavelength,
                                                     n2,
                                                     neff_max=10)
    setup.nonlinearity = [fun_gammax, fun_gammay]

    # The dispersion model is built from a OPD
    # read mat file for dispersion
    mat_path = os.path.join(os.path.dirname(__file__),
                            'data', 'laurent.mat')
    mat = gnlse.read_mat(mat_path)
    expL = mat['L'][0, 0]
    deltaOPD1 = mat['deltaOPD1'][:, 0]
    lambdasOPD1 = mat['lambdas1'][:, 0]
    OPD_lambda_um1 = mat['OPD_lambda_um1'][:, 0]
    deltaOPD2 = mat['deltaOPD2'][:, 0]
    lambdasOPD2 = mat['lambdas2'][:, 0]
    OPD_lambda_um2 = mat['OPD_lambda_um2'][:, 0]

    deltaN = 0.000336960439368341
    mat_path = os.path.join(os.path.dirname(__file__),
                            'data', 'loss_fited.mat')
    mat = gnlse.read_mat(mat_path)
    fiber_loss = mat['a'][0]  # dB/m
    lambdas_loss = mat['l'][0] * 1e3  # nm

    setup.dispersion_model = DubbleDispersionFiberFromOPD(
        lambdasOPD1, lambdasOPD2, deltaOPD1,
        deltaOPD2, expL, deltaN, w0, (lambdas_loss, fiber_loss))

    # Nonlinear Simulation
    ###########################################################################

    theta_list = [1, 89]
    count = len(theta_list)

    for i, th in enumerate(theta_list):
        theta = th / 180 * np.pi  # rad
        setup.pulse_model = DoubleSechEnvelope(
            input_energy, FWHM, theta, w0)

        solver = CNLSE(setup)
        solution = solver.run()

        # prepare initial vectors
        Z = solution.Z
        Atx = solution.At[:setup.z_saves, :]
        Aty = solution.At[setup.z_saves:, :]
        AWx = solution.AW[:setup.z_saves, :]
        AWy = solution.AW[setup.z_saves:, :]

        # filter noise & simulate futher distance
        if total_fiber_length > setup.fiber_length:
            lITx = 10 * np.log10(np.abs(Atx[-1, :]) ** 2,
                                 where=(np.abs(Atx[-1, :]) ** 2 > 0))
            lITy = 10 * np.log10(np.abs(Aty[-1, :]) ** 2,
                                 where=(np.abs(Aty[-1, :]) ** 2 > 0))
            lIT = lITx + lITy
            time_filter = solution.t[np.argmax(lIT)]
            iis = np.logical_or(
                solution.t < (time_filter - .5),
                solution.t > (time_filter + .5))
            Atx[-1, iis] = 0
            Aty[-1, iis] = 0

            # simulate futher distance
            for j in tqdm.tqdm(range(2, total_fiber_length + 1)):
                setup.pulse_model = np.concatenate(
                    (Atx[-1, :], Aty[-1, :]))
                solver = CNLSE(setup)
                solution = solver.run()
                Z = np.concatenate((Z, solution.Z[1:] + j - 1), axis=None)
                Atx = np.concatenate((Atx, solution.At[1:setup.z_saves, :]))
                Aty = np.concatenate((Aty, solution.At[setup.z_saves + 1:, :]))
                AWx = np.concatenate((AWx, solution.AW[1:setup.z_saves, :]))
                AWy = np.concatenate((AWy, solution.AW[setup.z_saves + 1:, :]))

            # update results
            solution.Z = Z
            solution.At = np.concatenate((Atx, Aty))
            solution.AW = np.concatenate((AWx, AWy))

            # save results
            solution.to_file(
                f'191119_polarisation_{total_fiber_length}m'
                f'_lambda_{setup.wavelength}nm'
                f'_power_{int(average_power * 1000)}mW'
                f'_angle_{int(theta / np.pi * 180)}deg_exp.mat')
