"""Dispersion operator in optical fibres.

Based on delivered frequency vector and damping indexes
script calculates linear dispersion operator in
frequency domain.

"""
import numpy as np
from scipy import interpolate

from gnlse.common import c
from gnlse.dispersion import Dispersion


class DubbleDispersionFiberFromOPD(Dispersion):
    """Calculates the dispersion in frequency domain

    Attributes
     ----------
    """

    def __init__(self, lambdasOPD1,
                 lambdasOPD2, deltaOPD1,
                 deltaOPD2, L, deltaN, w0,
                 loss=None, doping=0.18,
                 maxloss=1e2):
        self.lambdasOPD1 = lambdasOPD1
        self.lambdasOPD2 = lambdasOPD2
        self.deltaOPD1 = deltaOPD1
        self.deltaOPD2 = deltaOPD2
        self.L = L
        self.deltaN = deltaN
        self.w0 = w0
        if loss is None:
            self.doping = doping
            self.maxloss = maxloss
            self.lambda_loss = None
            self.fiber_loss = None
        elif loss == 0:
            self.doping = None
            self.maxloss = None
            self.lambda_loss = None
            self.fiber_loss = 0.
        else:
            self.doping = None
            self.maxloss = None
            self.lambda_loss = loss[0]
            self.fiber_loss = loss[1]

    def D(self, V):
        omegas1 = 2 * np.pi * c / (self.lambdasOPD1 * 1e3)  # rad*THz
        omegas2 = 2 * np.pi * c / (self.lambdasOPD2 * 1e3)  # rad*THz

        Neff1 = -self.deltaOPD1 * 1e-3 / self.L
        Neff2 = -self.deltaOPD2 * 1e-3 / self.L

        beta11 = Neff1 / c
        beta12 = Neff2 / c

        mu = np.mean(omegas1 - self.w0)
        std = np.std(omegas1 - self.w0, ddof=1)
        p11 = np.polyfit((omegas1 - self.w0 - mu)/std, beta11 * 1e9, 6)
        p11scaled = np.poly1d(p11)
        p1w01 = p11scaled(-mu/std)
        p01 = np.array(np.polyint(p11scaled))
        fp01 = np.poly1d(p01)
        p01[-1] = p01[-1] - fp01(-mu/std)
        beta01 = fp01((V - mu)/std) * std

        mu = np.mean(omegas2 - self.w0)
        std = np.std(omegas2 - self.w0, ddof=1)
        p12 = np.polyfit((omegas2 - self.w0 - mu)/std, beta12 * 1e9, 6)
        p12scaled = np.poly1d(p12)
        p1w02 = p12scaled(-mu/std)
        p02 = np.array(np.polyint(p12scaled))
        fp02 = np.poly1d(p02)
        p02[-1] = p02[-1] - fp02(-mu/std)
        beta02 = fp02((V - mu)/std) * std

        Bx = beta01 - (p1w01 + p1w02) / 2 * V
        By = beta02 - (p1w01 + p1w02) / 2 * V

        deltaB = self.deltaN * self.w0 / c * 1e9

        # Damping
        if self.fiber_loss is None:
            self.calc_loss(2 * np.pi * c / (V + self.w0) / 1e3)
        elif (
            isinstance(self.fiber_loss, int) or isinstance(self.fiber_loss, float)
            ) and self.fiber_loss == 0:
            self.alpha = 0
        else:
            WL = 2 * np.pi * c / (V + self.w0)  # nm
            # Extrapolate loss for a lambda vector
            loss_interp = interpolate.interp1d(self.lambda_loss,
                                               self.fiber_loss,
                                               kind='cubic',
                                               fill_value="extrapolate")
            loss = loss_interp(WL)
            loss[WL > self.lambda_loss[-1]] = np.max(self.fiber_loss)
            self.alpha = loss / (10 / np.log(10))

        # Linear dispersion operator
        Lx = 1j * Bx - self.alpha / 2
        Ly = 1j * By - self.alpha / 2

        return np.concatenate((Lx, Ly)), deltaB

    def calc_loss(self, lambdas):

        R = .74 + (2.33 - .74) * self.doping
        alphaR = R * lambdas**(-4) * 1e-3

        alphaoh_1_38 = 2.43

        sigma_lambda = 0.030
        alphaoh = 0.00012 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.444) / (sigma_lambda))**2) + \
                  0.00050 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.506) / (sigma_lambda))**2) + \
                  0.00030 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.566) / (sigma_lambda))**2) + \
                  0.00640 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.593) / (sigma_lambda))**2) + \
                  0.00028 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.651) / (sigma_lambda))**2) + \
                  0.00440 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.685) / (sigma_lambda))**2) + \
                  0.07800 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.724) / (sigma_lambda))**2) + \
                  0.00380 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.825) / (sigma_lambda))**2) + \
                  0.08000 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.878) / (sigma_lambda))**2) + \
                  1.6 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 0.943) / (sigma_lambda))**2) + \
                  0.07 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 1.139) / (sigma_lambda))**2) + \
                  2.7 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 1.246) / (sigma_lambda))**2) + \
                  alphaoh_1_38 * np.exp(-.5 * ((lambdas - 1.383) / (sigma_lambda))**2) + \
                  0.84 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 1.894) / (sigma_lambda))**2) + \
                  201 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 2.212) / (sigma_lambda))**2) + \
                  10000 / 62.7 * alphaoh_1_38 * np.exp(-.5 * ((lambdas - 2.722) / (sigma_lambda))**2)

        alphaIR = 4.2e8 * np.exp(-47.5 / lambdas)
        a = (alphaoh + alphaR + alphaIR) / (10 / np.log(10))
        a[a > self.maxloss] = self.maxloss
        self.alpha = a
