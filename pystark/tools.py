import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, k, c, atomic_mass, m_e
import pystark
import math

valid_x_units = ['m', 'Hz']


def get_approx_griem_a12_coefficient(n_upper):
    """
    Griem's scaling for density with FWHM is useful for simple analytical tests. Currently only has coefficients for H
    Balmer series.

    :param species:
    :param n_upper:
    :param n_lower:
    :return:
    """

    griem_a12s = {3: 0.05, 4: 0.08, 5: 0.0922747860122222, 6: 0.17, 7: 0.22, 8: 0.28, 9: 0.36,
                  10: 0.46}  # values approximate

    return griem_a12s[n_upper]


def get_wl_centre(species, n_upper, n_lower):
    """
    Central wavelength of spectral line ( m ). Source: NIST.

    Not fine-structure resolved (ie. weighted mean of fine structure components)

    :param species: string, 'H' etc.
    :param n_upper: int, upper principal quantum number
    :param n_lower: int, lower principal quantum number
    :return: float, wl_centre ( m )

    """

    wls_nm = {'H32': 656.279,
              'H42': 486.135,
              'H52': 434.0472,
              'H62': 410.1734,
              'H72': 397.0075,

              'D32': 656.106652,
              'D42': 486.00013,
              'D52': 433.92833,
              'D62': 410.06186,
              'D72': 396.89923,

              'He43': 468.6,
              'C00': 464.9,  # TODO fix this
           }
    key = species + str(n_upper) + str(n_lower)
    assert key in list(wls_nm.keys())
    return wls_nm[key] * 1e-9


def get_species_mass(species):
    """
    mass of emitting ion species ( kg )

    :param species: string, 'H', 'D', 'He' etc.
    :return: float, mass ( kg )
    """

    mass = {'H': 1.00794 * atomic_mass,
            'D': 2.01410178 * atomic_mass,
            'T': 3.01604928199 * atomic_mass,
            'He': 4.002602 * atomic_mass,
            'C': 12.0107 * atomic_mass,
            }
    assert species in list(mass.keys())
    return mass[species]


def get_wl_axis(species, n_upper, n_lower, dens, temp, bfield, no_fwhm=100, npts=100001, wl_centre=None):
    """
    For a given species, transition and plasma, return sensible, regular wavelength axis ( m )
    
    :param n_upper: 
    :param dens:
    :param temp:
    :param bfield:
    :param no_fwhm:
    :param npts:
    :param wl_centre:
    :return: 
    """

    if wl_centre is None:
        wl_centre = get_wl_centre(species, n_upper, n_lower)

    fwhm_voigt_hz = estimate_fwhm(species, n_upper, n_lower, dens, temp, bfield)
    fwhm_voigt_m = c * fwhm_voigt_hz / (c / wl_centre) ** 2
    return np.linspace(wl_centre - no_fwhm * fwhm_voigt_m, wl_centre + no_fwhm * fwhm_voigt_m, npts)


def estimate_fwhm(species, n_upper, n_lower, e_dens, temp, bfield):
    """
    use scalings to estimate the lineshape full width half maximum (FWHM)

    :param n_upper: upper principal quantum number
    :param e_dens: electron density [ /m^3 ]
    :param n_temp: neutral temperature [ eV ]
    :param bfield: magnetic field strength [ T ]
    :param isotope:

    :return: fwhm [ Hz ]

    """

    fwhm_doppler = get_fwhm_doppler(species, n_upper, n_lower, temp)
    fwhm_stark = get_fwhm_stark(species, n_upper, n_lower, e_dens)
    zeeman_split = e / (4 * np.pi * m_e) * bfield

    # total FWHM (Voigt approximation from https://en.wikipedia.org/wiki/Voigt_profile)
    fwhm = 0.5346 * fwhm_stark + np.sqrt(0.2166 * fwhm_stark ** 2 + fwhm_doppler ** 2)

    return fwhm + zeeman_split


def get_fwhm_doppler(species, n_upper, n_lower, temp):
    """
    Doppler-broadened FWHM
    
    :param n_upper: 
    :param temp: [ eV ]
    :param mass:

    :return: fwhm [ Hz ]

    """
    mass = get_species_mass(species)
    freq_centre = c / get_wl_centre(species, n_upper, n_lower)
    temp_k = temp * e / k
    v_th = np.sqrt(2 * k * temp_k / mass)
    sigma = v_th * freq_centre / (np.sqrt(2) * c)
    fwhm_hz = 2 * np.sqrt(2 * np.log(2)) * sigma

    return fwhm_hz


def get_fwhm_stark(species, n_upper, n_lower, e_dens):
    """
    Approximate Stark-broadened FWHM using tabulated values

    :param species:
    :param n_upper:
    :param n_lower:
    :param e_dens: [ m^-3 ] electron density

    :return: fwhm [ Hz ]
    """

    wl_centre = get_wl_centre(species, n_upper, n_lower)

    if species == 'He':
        """
        source: He II I' Stark broadening and intensity ratio of C IV and C III lines calibrated
        with Thomson scattering for high-density plasma diagnostics, A. Gawran et al., phys rev. A (1988)
        """

        fwhm_nm = 3.65e-20 * e_dens ** 0.826
        fwhm_hz = fwhm_nm * 1e-9 * c / wl_centre ** 2

    else:
        alpha_12 = get_approx_griem_a12_coefficient(n_upper)
        ne_20 = e_dens * 1e-20  # convert into units: 10 ** 20 m ** -3
        # fwhm_m = 1e-9 * 0.54 * alpha_12 * ne_20 ** (2 / 3)  # Griem scaling
        fwhm_m = 1e-9 * 0.53860867250797 * alpha_12 * ne_20 ** (2 / 3)  # Griem scaling
        fwhm_hz = fwhm_m * c / wl_centre ** 2

    return fwhm_hz


def get_stehle_balmer_wavelength(n_upper):
    # ensure given n_upper + n_lower fall within tabulated values
    n_lower = 2
    assert n_upper in range(n_lower + 1, 31)

    prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'

    # extract raw tabulated tabulated_data
    # tab_temp_k = np.array(pystark.nc.variables[prefix + 'tempe'].tabulated_data)  # tabulated electron temperatures (K)
    olam0, = pystark.nc.variables[prefix + 'olam0'].data  # line centre wavelength (A)
    olam0 *= 1e-10

    return olam0


def get_fwhm(x, y, disp=False):
    """
    given a function with a SINGLE PEAK, find the FWHM without fitting. QUICK AND DIRTY.
    """

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


def ls_norm(ls, x, norm_type='area'):
    """ lineshape normalisation. """

    valid_norm_types = ['area', 'height']
    assert norm_type in valid_norm_types

    if norm_type == 'area':
        return ls / np.trapz(ls, x)

    if norm_type == 'height':
        return ls / np.max(ls)


def FINTRP(x1, x2, x3, y1, y2, y3, x):
    """
    Interpolation routine used in Stehle calculation.
    """

    if x == x2:
        return y2

    a12 = x1 - x2
    a22 = x1 - x3
    v1 = y1 - y2
    v2 = y1 - y3

    if ((y1 < y2) and (y2 < y3)) or ((y1 > y2) and (y2 > y3)):
        deter = v1 * a22 - v2 * a12
        if np.abs(deter) < 1.0E-40:
            return y1 + (x - x1) * (y3 - y1) / (x3 - x1)
        a21 = x1 * y1
        a11 = a21 - x2 * y2
        a21 = a21 - x3 * y3
        c = (a22 * a11 - a12 * a21) / deter
        a = (-v2 * a11 + v1 * a21) / deter
        b = (y1 - a) * (x1 - c)
        return a + b / (x - c)
    else:
        x1c = x1 * x1
        a11 = x1c - x2 * x2
        a21 = x1c - x3 * x3
        deter = a11 * a22 - a12 * a21
        if np.abs(deter) < 1.0E-40:
            raise Exception('FINTRP error: incorrect inputs')
        a = (a22 * v1 - a12 * v2) / deter
        b = (-a21 * v1 + a11 * v2) / deter
        c = y1 - a * x1c - b * x1
        return (a * x + b) * x + c


def to_precision(x, p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0. and p == 0:
        return "0"

    if x == 0.:
        return "0." + "0" * (p - 1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x / tens)

    if n < math.pow(10, p - 1):
        e = e - 1
        tens = math.pow(10, e - p + 1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens - x):
        n = n + 1

    if n >= math.pow(10, p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p - 1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e + 1])
        if e + 1 < len(m):
            out.append(".")
            out.extend(m[e + 1:])
    else:
        out.append("0.")
        out.extend(["0"] * -(e + 1))
        out.append(m)

    return "".join(out)


def doppler_lineshape(species, x, x_centre, temp, x_units='m'):
    """
    generate Doppler broadened lineshape, area-normalised to 1.
    
    Assumed Maxwellian velocity distribution -- Gaussian profile.
    
    :param wavelengths: [ m ]
    :param wavelength_centre: [ m ]
    :param temp: [ eV ]
    :return: 
    """
    mass = get_species_mass(species)
    temp_k = temp * e / k  # temperature [ K ]
    v_th = np.sqrt(2 * k * temp_k / mass)  # thermal speed [ m/s ]

    assert x_units in valid_x_units

    if x_units == 'm':
        freqs, freq_centre = convert_ls_units(x, x_centre, mode='direct')
    else:
        freqs, freq_centre = x, x_centre

    sigma = v_th * freq_centre / (np.sqrt(2) * c)  # [ Hz ]
    ls_d = (freq_centre ** -1) * np.sqrt((c / v_th) ** 2 / np.pi) * np.exp(- 0.5 * ((freqs - freq_centre) / sigma) ** 2)

    return ls_d


def get_freq_axis(species, n_upper, n_lower, e_dens, temp, bfield, no_fwhm=50, npts=5001, wl_centre=None):
    """

    :param n_upper:
    :param e_dens:
    :param n_temp:
    :param bfield:
    :param no_fwhm:
    :param npts:
    :param wl_centre:
    :return:
    """

    if wl_centre is None:
        wl_centre = get_wl_centre(species, n_upper, n_lower)

    freq_centre = c / wl_centre
    fwhm_hz = estimate_fwhm(species, n_upper, n_lower, e_dens, temp, bfield)

    return np.linspace(freq_centre - no_fwhm * fwhm_hz, freq_centre + no_fwhm * fwhm_hz, npts)


def get_freq_axis_conv(freq_axis, extra=1000):
    """
    extend a uniform frequency axis for use in convolution to avoid edge effects

    """

    min_freq, max_freq, = np.min(freq_axis), np.max(freq_axis)
    len_axis = len(freq_axis)
    dfreq = abs(freq_axis[1] - freq_axis[0])
    min_freq_conv = min_freq - extra / 2 * dfreq
    max_freq_conv = max_freq + extra / 2 * dfreq
    freq_axis_conv = np.linspace(min_freq_conv, max_freq_conv, len_axis + extra)

    return freq_axis_conv


def convert_ls_units(x, x_centre, mode='uniform', x_out=None, ls=None):
    """ convert a lineshape between frequency / wavelength spaces.
    
    output axis is in ascending order (so has been flipped)
    
    :param x: 
    :param x_centre: 
    :param mode: does the output grid need to be uniform? 
    :param ls: optional, convert the lineshape values themselves
    :return: 
    """

    # preliminary checks
    valid_modes = ['direct', 'uniform', 'interp']
    assert mode in valid_modes
    assert np.min(x) < x_centre < np.max(x)

    # if an output axis is supplied, activate interp mode automatically
    if x_out is not None:
        mode = 'interp'
    # if interp mode selected and a lineshape not supplied, raise an error
    if mode == 'interp' and ls is None:
        raise Exception()

    x_centre_out = c / x_centre

    # convert x grid values
    if mode == 'direct':
        x_out = c / x[::-1]

    elif mode == 'uniform':
        x_min, x_max = np.min(x), np.max(x)
        npts = len(x)
        x_out = np.linspace(c / x_max, c / x_min, npts)

    elif mode == 'interp':
        # need to interpolate ls intensities onto reciprocal of output grid
        ls = np.interp(c / x_out[::-1], x, ls)

    # convert lineshape intensities
    if ls is None:
        return x_out, x_centre_out
    else:
        # assert len(x) == len(ls)
        ls_out = ls[::-1] * c * x_out ** -2  # conservation of energy!
        return x_out, x_centre_out, ls_out


def find_nearest_idx(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()