import numpy as np
import matplotlib.pyplot as plt
import pystark
from scipy.constants import *
from scipy.signal import fftconvolve as conv


def simple_profile(n_upper, n_lower, temperature, density, wavelengths, model='lomanowski', display=False):
    """ Approximate Stark-broadened line profile, area-normalised to 1.

     use mode='griem' for a Lorentzian Stark profile that follows Griem's density scaling:

    - Scaling is well fulfilled for n_e > 10 ** 20 m ** -3 and T_i ~ 1eV.
    - Unverified for n_e > 10 ** 21 m ** -3


    use mode='lomanowski' for B. Lomanowski's parameterised Stark profile, based on Stehle's tabulated tabulated_data.

    :param n_upper: upper principal quantum number
    :param n_lower: lower principal quantum number
    :param temperature: in eV
    :param density: in m ** -3
    :param spectrum: output the spectrum as a function of 'frequency' or 'wavelength'.

    """

    valid_n_upper = [3, 4, 5, 6, 7, 8, 10]
    assert n_upper in valid_n_upper

    valid_modes = ['griem', 'lomanowski']
    assert model in valid_modes

    # line centre wavelength from NIST
    # prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'
    # wl_0 = pystark.nc.variables[prefix + 'olam0'].tabulated_data[0] * 1e-10  # line centre wavlength (m)
    wl_0 = pystark.get_NIST_balmer_wavelength(n_upper)
    wl_0_nm = wl_0 * 1e9
    freq_0 = c / wl_0

    # ----- approximate FWHM of Voigt profile to dynamically generate frequency axis
    fwhm_voigt_hz, fwhm_lorentz_hz, fwhm_gauss_hz = pystark.estimate_fwhms(n_upper, density, temperature)
    sigma_gauss_hz = fwhm_gauss_hz / (2 * np.sqrt(2 * np.log(2)))

    # generate frequency and wavelength axes
    num_fwhm = 20
    len_axis = 2001  # odd number of points such that nu_0 lies exactly at the array centre.
    min_freq, max_freq = (freq_0 - num_fwhm * fwhm_voigt_hz, freq_0 + num_fwhm * fwhm_voigt_hz)
    freq_axis = np.linspace(min_freq, max_freq, len_axis)

    min_wl, max_wl = (c / max_freq, c / min_freq)
    wl_axis = np.linspace(min_wl, max_wl, len_axis)

    # generate a larger, convolution axis -- this avoids any
    # unwanted edge effects from the convolution.

    dfreq = (max_freq - min_freq) / len_axis
    min_freq_conv = min_freq - len_axis / 3 * dfreq
    max_freq_conv = max_freq + len_axis / 3 * dfreq
    freq_axis_conv = np.linspace(min_freq_conv, max_freq_conv, 3 * len_axis)

    dwl = (max_wl - min_wl) / len_axis
    min_wl_conv = min_wl - len_axis / 3 * dwl
    max_wl_conv = max_wl + len_axis / 3 * dwl
    wl_axis_conv = np.linspace(min_wl_conv, max_wl_conv, 3 * len_axis)
    wl_axis_conv_nm = wl_axis_conv * 1e9

    # Doppler lineshape
    temperature_k = temperature * e / k
    mass_kg = physical_constants['deuteron mass'][0]
    v_th = np.sqrt(2 * k * temperature_k / mass_kg)
    doppler_lineshape_hz = (freq_0 ** -1) * ((np.pi * (v_th ** 2 / (c ** 2))) ** -0.5) * np.exp \
        (-(((freq_axis_conv - freq_0) ** 2) / (2 * sigma_gauss_hz ** 2)))

    # Stark lineshape
    if model == 'griem':

        hwhm_lorentz_hz = fwhm_lorentz_hz / 2
        stark_lineshape_hz = (1 / np.pi) * hwhm_lorentz_hz / ((freq_axis_conv - freq_0) ** 2 + hwhm_lorentz_hz ** 2)

    elif model == 'lomanowski':

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

        s_lineshape_m = 1 / (abs(wl_axis_conv_nm - wl_0_nm) ** (5. / 2.) + (delta_lambda_12ij / 2) ** (5. / 2.))
        s_lineshape_m /= np.trapz(s_lineshape_m, wl_axis_conv)

        stark_lineshape_hz = np.interp(freq_axis_conv, c / wl_axis_conv[::-1], s_lineshape_m[::-1] * (wl_axis_conv ** 2 / c))

    # Voigt lineshape
    # voigt_lineshape_hz = np.convolve(stark_lineshape_hz, doppler_lineshape_hz, 'same')
    voigt_lineshape_hz = conv(stark_lineshape_hz, doppler_lineshape_hz, 'same')
    voigt_lineshape_hz /= np.trapz(voigt_lineshape_hz, freq_axis_conv)  # normalise

    d_lineshape_m = np.interp(wavelengths, c / freq_axis_conv[::-1], doppler_lineshape_hz[::-1] * freq_axis_conv ** 2 / c)
    s_lineshape_m = np.interp(wavelengths, c / freq_axis_conv[::-1], stark_lineshape_hz[::-1] * freq_axis_conv ** 2 / c)
    sd_lineshape_m = np.interp(wavelengths, c / freq_axis_conv[::-1], voigt_lineshape_hz[::-1] * freq_axis_conv ** 2 / c)

    if display:

        # plot
        fsize=14
        wavelengths_nm = wavelengths * 1e9
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.plot(wavelengths_nm, d_lineshape_m, label='Doppler')
        ax1.plot(wavelengths_nm, s_lineshape_m, label='Stark')
        ax1.plot(wavelengths_nm, sd_lineshape_m, label='Stark-Doppler')

        leg = ax1.legend(fontsize=fsize)
        ax1.set_xlabel('wavelength (nm)', size=fsize)
        ax1.set_ylabel('normalised lineshape', size=fsize)
        plt.show()

    return sd_lineshape_m, s_lineshape_m, d_lineshape_m

