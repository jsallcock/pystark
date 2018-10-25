import pystark
import numpy as np
from scipy.constants import c
import matplotlib.pyplot as plt


def demo():
    """ Demo script compares all available lineshape models for selected Balmer line. """

    n_upper = 6
    n_lower = 2
    temp = 5.  # [eV]
    dens = 5e20  # [m-3]
    bfield = 1.5  # [T]
    viewangle = 0  # [deg]

    # generate appropriate wavelength axis

    wl_0 = pystark.get_NIST_balmer_wavelength(n_upper)
    fwhm_voigt_hz, _, _ = pystark.estimate_fwhms(n_upper, dens, temp)
    fwhm_voigt_m = c * fwhm_voigt_hz / (c / wl_0) ** 2

    wl_axis = pystark.generate_wavelength_axis(n_upper, temp, dens)
    wl_axis_nm = wl_axis * 1e9

    display=False

    # generate plots
    fig = plt.figure(figsize=[6, 6])
    fsize=14
    ax = fig.add_subplot(111)
    ax.set_title('$n_e = $' + str(dens) + ' m$^{-3}$ \n $T = $' + str(temp) + ' eV \n $B = $' + str(bfield) + ' T')

    try:
        ls_stehle = pystark.stehle_profile(n_upper, n_lower, temp, dens, wl_axis, display=display)
        ax.plot(wl_axis_nm, pystark.norm_a(ls_stehle, wl_axis), '-', label='Stehle')
    except Exception as e:
        print('--Stehle calculation failed:')
        print(e)

    try:
        ls_rosato = pystark.rosato_profile(n_upper, dens, temp, bfield, viewangle, wl_axis, display=display)
        ax.plot(wl_axis_nm, pystark.norm_a(ls_rosato, wl_axis), '-', label='Rosato')
    except Exception as e:
        print('--Rosato calculation failed:')
        print(e)

    try:
        ls_loman, _, _ = pystark.simple_profile(n_upper, n_lower, temp, dens, wl_axis, model='lomanowski')
        ax.plot(wl_axis_nm, pystark.norm_a(ls_loman, wl_axis), '-', label='Lomanowski')
    except Exception as e:
        print('--Lomanowski calculation failed:')
        print(e)

    try:
        ls_griem, _, _ = pystark.simple_profile(n_upper, n_lower, temp, dens, wl_axis, model='griem')
        ax.plot(wl_axis_nm, pystark.norm_a(ls_griem, wl_axis), '-', label='Griem')
    except Exception as e:
        print('--Griem calculation failed:')
        print(e)

    ax.set_xlim([np.min(wl_axis_nm), np.max(wl_axis_nm)])
    leg = ax.legend(fontsize=fsize)
    ax.set_xlabel('wavelength (nm)', size=fsize)
    plt.show()

    return

if __name__ == '__main__':
    demo()
