import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, c, k, h, physical_constants
from scipy.signal import fftconvolve as conv
from pystark import rosato_path, rosato_database_path
import pystark

sys.path.insert(0, rosato_path)
import LS_DATA_read_f2py  # this is the shared object generated using the fortran subroutines in LS_DATA_read_f2py.f90


def rosato_profile(n_upper, dens, temp, bfield, viewangle, wavelengths, line_centre=None, display=False):
    """ Stark-Zeeman-Doppler lineshape profile for the first Balmer lines, interpolated from the Rosato et al. tabulated_data 
    tables.
    
    Ported by Joe Allcock, 08/18.
    
    :param n_upper: upper principal quantum number, must be in [3, 7]
    :param dens: [m-3] density, must be in [10^19, 10^22]
    :param temp: [eV] temperature, must be in [0.316, 31.6]
    :param bfield: [T] magnetic field strength, must be in [0, 5[
    :param viewangle: [deg]
    :param wavelengths: [m] input wavelength axis
    :param line_centre: [m] line centre wavelength, optional.
    :param display: bool.
    
    :return: 
    """

    dens *= 1e-6

    if line_centre is None:
        # Use the NIST value for chosen Balmer line.
        centre_m = pystark.get_NIST_balmer_wavelength(n_upper)
    else:
        centre_m = line_centre

    centre_ev = c * h / (e * centre_m)
    centre_hz = c / centre_m

    # wavelengths [m] -> [eV]
    interp_wavelengths_ev = c * h / (e * wavelengths)

    wmax = np.max(abs(interp_wavelengths_ev - centre_ev))  # [ev] furthest point from the centre required in the interpolation
    npts = 1001  # must be odd

    ##### Stark-Zeeman lineshape #####

    detunings_sz_ev, ls_sz_ev = rosato_stark_zeeman_profile(n_upper, dens, temp, bfield, viewangle, wmax, npts)
    detunings_sz_hz = e * detunings_sz_ev / h
    ls_sz_hz = ls_sz_ev * h / e

    # generate a larger, convolution axis -- avoiding any
    # unwanted edge effects from the convolution.

    d_ev = detunings_sz_ev[-1] - detunings_sz_ev[-2]
    num_extra_points = 100
    detunings_sz_ev_extended = np.linspace((-wmax - (num_extra_points / 2) * d_ev),
                                           (wmax + (num_extra_points / 2) * d_ev), npts + num_extra_points)
    detunings_sz_hz_extended = e * detunings_sz_ev_extended / h

    ##### Doppler lineshape #####
    temp_k = temp * e / k
    mass_kg = physical_constants['deuteron mass'][0]
    v_th = np.sqrt(2 * k * temp_k / mass_kg)
    sigma_gauss_hz = v_th * centre_hz / (np.sqrt(2) * c)
    ls_d_hz = (centre_hz ** -1) * ((np.pi * (v_th ** 2 / (c ** 2))) ** -0.5) * np.exp(-(((detunings_sz_hz_extended) ** 2) / (2 * sigma_gauss_hz ** 2)))

    # convolution -- in frequency space
    ls_szd_hz = conv(ls_sz_hz, ls_d_hz, 'same')
    detunings_szd_hz = detunings_sz_hz
    ls_szd_hz /= np.trapz(ls_szd_hz, detunings_szd_hz)  # normalise

    # convert to wavelength
    detunings_szd_m = c * detunings_szd_hz / (detunings_szd_hz + centre_hz) ** 2
    ls_szd_m = c * ls_szd_hz / (detunings_szd_m + centre_m) ** 2

    # interpolate onto input axis
    ls_szd_m_interp = np.interp(wavelengths, detunings_szd_m + centre_m, ls_szd_m)

    if display:
        # interpolate other components onto input axis for display

        # Stark-Zeeman
        detunings_sz_m = c * detunings_sz_hz / (detunings_sz_hz + centre_hz) ** 2
        ls_sz_m = c * ls_sz_hz / (detunings_sz_m + centre_m) ** 2
        ls_sz_m_interp = np.interp(wavelengths, detunings_sz_m + centre_m, ls_sz_m)

        # Doppler
        detunings_sz_m_extended = c * detunings_sz_hz_extended / (detunings_sz_hz_extended + centre_hz) ** 2
        ls_d_m = c * ls_d_hz / (detunings_sz_m_extended + centre_m) ** 2
        ls_d_m_interp = np.interp(wavelengths, detunings_sz_m_extended + centre_m, ls_d_m)

        # plot
        fsize=14
        wavelengths_nm = wavelengths * 1e9
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.plot(wavelengths_nm, ls_szd_m_interp, label='Stark-Zeeman-Doppler')
        ax1.plot(wavelengths_nm, ls_sz_m_interp, label='Stark-Zeeman')
        ax1.plot(wavelengths_nm, ls_d_m_interp, label='Doppler')

        leg = ax1.legend(fontsize=fsize)
        ax1.set_xlabel('wavelength (nm)', size=fsize)
        ax1.set_ylabel('normalised lineshape', size=fsize)
        plt.show()


    return ls_szd_m_interp


def rosato_stark_zeeman_profile(n_upper, dens, temp, bfield, viewangle, wmax, npts, display=False):
    """
     Stark-Zeeman lineshape profile for the first Balmer lines, interpolated from the Rosato et al. tabulated_data.
    
    Essentially a python wrapper for the LS_READ_data.f90
    
    :param n_upper: upper principal quantum number
    :param dens: [cm-3]
    :param temp: [eV]
    :param bfield: [T]
    :param viewangle: [deg] 
    :param wmax: [eV]
    :param npts: 
    :param display: bool.
    :return: 
    """
    balmer_line_names = ['-', '-', '-', 'D_alpha', 'D_beta', 'D_gamma', 'D_delta', 'D_epsilon']
    # list index gives transition n_upper

    detunings_axis = np.linspace(-wmax, wmax, npts)

    # overwrite the in.txt fortran input file -- is this necessary?
    # dir_in = os.path.join(rosato_path).encode('UTF-8')
    # dir_in = str(os.path.join(rosato_path))
    # dir_in = 'hello world'.encode('utf-8')
    # dir_in = 'hello world from python'
    # # print(len(dir_in), dir_in)


    # get the parameter grid bound indices
    iN, iT, iB = LS_DATA_read_f2py.set_bounds(dens, temp, bfield)

    viewangle_idxs = [0, 1]
    viewangle_rad = viewangle * np.pi / 180
    ls = np.zeros([npts, 2])

    for i, viewangle_idx in enumerate(viewangle_idxs):
        name = LS_DATA_read_f2py.set_name_file(n_upper, iN, iT, iB, viewangle_idx)

        # replace the original 'dir' outputted from LS_DATA_read with one that is defined using the full path
        # -- this step could cause headaches later on, but allows fn to be called from anywhere.

        dir = os.path.join(rosato_database_path, balmer_line_names[n_upper] + '/').encode('utf-8')

        w_arr, ls_arr = LS_DATA_read_f2py.read_file(len(dir), dir, name)
        ls[:, i] = LS_DATA_read_f2py.ls_interpol(n_upper, dens, temp, bfield, wmax, npts, w_arr, ls_arr, iN, iT, iB)

    ls = ls[:, 0] * np.sin(viewangle_rad) ** 2 + ls[:, 1] * np.cos(viewangle_rad) ** 2

    if display:
        plt.figure()
        plt.plot(detunings_axis, ls)
        plt.ylabel('ls')
        plt.xlabel('Detuning (eV)')
        # plt.semilogy()
        plt.show()

    return detunings_axis, ls



if __name__ == '__main__':
    s = time.time()
    wls = np.linspace(656, 657, 1000) * 1e-9
    rosato_profile(3, 5.e19, 3.16, 4.99, 90., wls, display=False)
    e = time.time()
    # print(e - s, 'sec')
    # x, y = rosato_stark_zeeman_profile(3, 2.41e15, 3., 5., 86.2, 1.e-1, 10001, display=True)