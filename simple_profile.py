import scipy.io.netcdf as netcdf
from scipy.constants import c, e, k
import matplotlib.pyplot as plt
import numpy as np
import time
import pystarky
from scipy.constants import *
import os
import scipy.integrate
from matplotlib.colors import LogNorm


def simple_profile(n_upper, n_lower, temperature, density, spectrum='wavelength', mode='griem'):
    """ Use Griem's scaling for an approximate Stark-broadened line in a certain parameter region. 

    - Scaling is well fulfilled for n_e > 10 ** 20 m ** -3 and T_i ~ 1eV.
    - Unverified for n_e > 10 ** 21 m ** -3

    TODO: CORRECT NORMALISATION TO MATCH STEHLE

    :param n_upper: upper principal quantum number
    :param n_lower: lower principal quantum number
    :param temperature: in eV
    :param density: in m ** -3

    """

    valid_n_upper = [4, 5, 6, 8, 10]

    # lomanowski_valid_n_upper = []
    assert n_upper in valid_n_upper

    valid_modes = ['griem', 'lomanowski']
    assert mode in valid_modes

    # load line centre from Stehle tables
    prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'
    wl_0 = nc.variables[prefix + 'olam0'].data[0] * 1e-10  # line centre wavlength (m)
    wl_0_nm = wl_0 * 1e9
    freq_0 = c / wl_0

    # ----- approximate FWHM of Voigt profile to dynamically generate frequency axis

    # Doppler FWHM: Gaussian
    temperature_k = temperature * e / k
    mass_kg = physical_constants['deuteron mass'][0]
    v_th = np.sqrt(2 * k * temperature_k / mass_kg)
    sigma_gauss_hz = v_th * freq_0 / (np.sqrt(2) * c)
    fwhm_gauss_hz = 2 * np.sqrt(2 * np.log(2)) * sigma_gauss_hz

    # Stark FWHM: Lorentzian
    griem_alpha_12s = {4: 0.08, 5: 0.09, 6: 0.17, 8: 0.28, 10: 0.46}
    alpha_12 = griem_alpha_12s[n_upper]

    ne_20 = density * 1e-20  # convert into units: 10 ** 20 m ** -3
    fwhm_lorentz_m = 1e-9 * 0.54 * alpha_12 * ne_20 ** (2. / 3.)  # Griem scaling
    fwhm_lorentz_hz = fwhm_lorentz_m * freq_0 ** 2 / c

    # total FWHM: Voigt
    fwhm_voigt_hz = 0.5346 * fwhm_lorentz_hz + np.sqrt(0.2166 * fwhm_lorentz_hz ** 2 + fwhm_gauss_hz ** 2)

    # generate frequency and wavelength axes
    num_fwhm = 20
    num_points = 20001  # odd number of points such that nu_0 lies exactly at the array centre.
    min_freq, max_freq = (freq_0 - num_fwhm * fwhm_voigt_hz, freq_0 + num_fwhm * fwhm_voigt_hz)
    freq_axis = np.linspace(min_freq, max_freq, num_points)

    min_wl, max_wl = (c / max_freq, c / min_freq)
    wl_axis = np.linspace(min_wl, max_wl, num_points)
    wl_axis_nm = wl_axis * 1e9

    # Doppler lineshape
    doppler_lineshape_hz = (freq_0 ** -1) * ((np.pi * (v_th ** 2 / (c ** 2))) ** -0.5) * np.exp \
        (-(((freq_axis - freq_0) ** 2) / (2 * sigma_gauss_hz ** 2)))

    # Stark lineshape
    if mode == 'griem':

        hwhm_lorentz_hz = fwhm_lorentz_hz / 2
        stark_lineshape_hz = (1 / np.pi) * hwhm_lorentz_hz / ((freq_axis - freq_0) ** 2 + hwhm_lorentz_hz ** 2)

    elif mode == 'lomanowski':

        # Paramaterised MMM Stark profile coefficients from Bart's paper
        loman_abc_ij_idx = {'32': 0,
                            '42': 1,
                            '52': 2,
                            '62': 3,
                            '72': 4,
                            '82': 5,
                            '92': 6,
                            '43': 7,
                            '53': 8,
                            '63': 9,
                            '73': 10,
                            '83': 11,
                            '93': 12}
        loman_c_ij = [3.710e-18,
                      8.425e-18,
                      1.310e-15,
                      3.954e-16,
                      6.258e-16,
                      7.378e-16,
                      8.947e-16,
                      1.330e-16,
                      6.640e-16,
                      2.481e-15,
                      3.270e-15,
                      4.343e-15,
                      5.588e-15]

        loman_a_ij = [0.7665,
                      0.7803,
                      0.6796,
                      0.7149,
                      0.7120,
                      0.7159,
                      0.7177,
                      0.7449,
                      0.7356,
                      0.7118,
                      0.7137,
                      0.7133,
                      0.7165]

        loman_b_ij = [0.064,
                      0.050,
                      0.030,
                      0.028,
                      0.029,
                      0.032,
                      0.033,
                      0.045,
                      0.044,
                      0.016,
                      0.029,
                      0.032,
                      0.033]

        ij_idx = loman_abc_ij_idx[str(n_upper) + str(n_lower)]
        c_ij = loman_c_ij[ij_idx]
        a_ij = loman_a_ij[ij_idx]
        b_ij = loman_b_ij[ij_idx]

        delta_lambda_12ij = c_ij * (density ** a_ij) / (temperature ** b_ij)  # nm

        stark_lineshape_m = 1 / (abs(wl_axis_nm - wl_0_nm) ** (5. / 2.) + (delta_lambda_12ij / 2) ** (5. / 2.))
        stark_lineshape_m /= np.trapz(stark_lineshape_m, wl_axis)

        stark_lineshape_hz = np.interp(freq_axis, c / wl_axis[::-1], stark_lineshape_m[::-1] * (wl_axis ** 2 / c))

    # Voigt lineshape
    voigt_lineshape_hz = np.convolve(doppler_lineshape_hz, stark_lineshape_hz, 'same')
    voigt_lineshape_hz /= np.trapz(voigt_lineshape_hz, freq_axis)  # normalise


    # return spectra in wavelength or frequency:
    assert spectrum in ['wavelength', 'frequency']

    if spectrum == 'wavelength':

        doppler_lineshape_m = np.interp(wl_axis, c / freq_axis[::-1], doppler_lineshape_hz[::-1] * freq_axis ** 2 / c)
        stark_lineshape_m = np.interp(wl_axis, c / freq_axis[::-1], stark_lineshape_hz[::-1] * freq_axis ** 2 / c)
        voigt_lineshape_m = np.interp(wl_axis, c / freq_axis[::-1], voigt_lineshape_hz[::-1] * freq_axis ** 2 / c)

        return wl_axis, wl_0, voigt_lineshape_m, stark_lineshape_m, doppler_lineshape_m

    elif spectrum == 'frequency':
        return freq_axis, freq_0, voigt_lineshape_hz, stark_lineshape_hz, doppler_lineshape_hz
