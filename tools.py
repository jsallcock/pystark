import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
import pystark

def find_nearest_idx(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def generate_frequency_axis(n_upper, n_lower, temperature, density):
    """ For a given transition and plasma parameters, return a regular frequency axis with sensible bounds using a voigt
     approximation for the line's FWHM. """

    freq_0 = c / get_NIST_balmer_wavelength(n_upper)
    fwhm_voigt_hz, fwhm_lorentz_hz, fwhm_gauss_hz = estimate_fwhms(n_upper, density, temperature)
    num_fwhm = 10
    num_points = 5001  # odd number of points such that nu_0 lies exactly at the array centre.
    min_freq, max_freq = (freq_0 - num_fwhm * fwhm_voigt_hz, freq_0 + num_fwhm * fwhm_voigt_hz)
    frequency_axis = np.linspace(min_freq, max_freq, num_points)


    return frequency_axis


def estimate_fwhms(n_upper, density, temperature):
    """ Estimate the FWHM for Stark and Doppler broadened Balmer line -- used for generating dynamic frequency detuning
    axes. """
    freq_0 = c / get_NIST_balmer_wavelength(n_upper)

    # Doppler FWHM: Gaussian
    temperature_k = temperature * e / k
    mass_kg = physical_constants['deuteron mass'][0]
    v_th = np.sqrt(2 * k * temperature_k / mass_kg)
    sigma_gauss_hz = v_th * freq_0 / (np.sqrt(2) * c)
    fwhm_gauss_hz = 2 * np.sqrt(2 * np.log(2)) * sigma_gauss_hz

    # Stark FWHM: Lorentzian
    griem_alpha_12s = {3:0.05, 4: 0.08, 5: 0.09, 6: 0.17, 7:0.22, 8: 0.28, 9:0.36, 10: 0.46}  # odd values are bogus, even values from Hutchinson -- this is just approximate though
    alpha_12 = griem_alpha_12s[n_upper]

    ne_20 = density * 1e-20  # convert into units: 10 ** 20 m ** -3
    fwhm_lorentz_m = 1e-9 * 0.54 * alpha_12 * ne_20 ** (2. / 3.)  # Griem scaling
    fwhm_lorentz_hz = fwhm_lorentz_m * freq_0 ** 2 / c

    # total FWHM: Voigt
    fwhm_voigt_hz = 0.5346 * fwhm_lorentz_hz + np.sqrt(0.2166 * fwhm_lorentz_hz ** 2 + fwhm_gauss_hz ** 2)

    return fwhm_voigt_hz, fwhm_lorentz_hz, fwhm_gauss_hz


def get_NIST_balmer_wavelength(n_upper):
    """ Source: NIST. These values are for deuterium, with unresolved fine-structure transitions/"""
    assert n_upper in [3, 4, 5, 6, 7, 8]
    return [0, 0, 0, 656.1012, 486.00013, 433.92833, 410.06186, 396.89923, 388.79902][n_upper] * 1e-9

def get_stehle_balmer_wavelength(n_upper):
    # ensure given n_upper + n_lower fall within tabulated values
    n_lower = 2
    assert n_upper in range(n_lower + 1, 31)

    prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'

    # extract raw tabulated data
    # tab_temp_k = np.array(pystark.nc.variables[prefix + 'tempe'].data)  # tabulated electron temperatures (K)
    olam0, = pystark.nc.variables[prefix + 'olam0'].data  # line centre wavelength (A)

    olam0 *= 1e-10

    return olam0


def get_fwhm(x, y, disp=False):
    """ given a function with a SINGLE PEAK, find the FWHM without fitting. QUICK AND DIRTY. """

    # normalise height
    y_norm = y / np.max(y)

    # split array into two about the max value
    idx_max = np.argmax(y_norm)
    x_l, x_u, y_l, y_u = x[:idx_max], x[idx_max+1:], y_norm[:idx_max], y_norm[idx_max+1:]

    # 1D interpolation
    hm_idx_l, hm_idx_u = np.interp(0.5, y_l, x_l), np.interp(0.5, y_u[::-1], x_u[::-1])

    fwhm = hm_idx_u - hm_idx_l


    if disp:

        print(fwhm)
        plt.figure()
        plt.plot(x, y_norm, '-')
        plt.plot(hm_idx_l, 0.5, '.')
        plt.plot(hm_idx_u, 0.5, '.')

        plt.show()

    return fwhm

def norm_h(ls):
    """ height normalise peak to 1. """
    return ls / np.max(ls)

def norm_a(ls, x):
    """ area normalise peak to 1. """
    return ls / np.trapz(ls, x)



if __name__ == '__main__':
    get_stehle_balmer_wavelength(5)