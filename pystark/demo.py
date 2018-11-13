import pystark
from pystark.tools import to_precision as tp

import numpy as np
import os, sys, traceback, time

from scipy.constants import c
import matplotlib.pyplot as plt


def demo():
    """ Demo script compares all available lineshape models for selected Balmer line. """

    # specify plasma
    temp = 5.  # [eV]
    dens = 2e20  # [m-3]
    bfield = 4.9  # [T]
    viewangle = 0 * np.pi / 180  # [rad]

    # specify line
    isotope = 'D'
    n_upper = 6
    n_lower = 2

    # plot params
    norm_type = 'area'
    fsize = 14

    # generate appropriate wavelength axis

    wl_0 = pystark.get_wl_centre(n_upper)
    fwhm_voigt_hz = pystark.estimate_fwhm(n_upper, dens, temp, bfield, isotope)
    fwhm_voigt_m = c * fwhm_voigt_hz / (c / wl_0) ** 2

    wl_axis = pystark.get_wavelength_axis(n_upper, dens, temp, bfield)
    wl_axis_nm = wl_axis * 1e9

    # generate plots
    fig = plt.figure(figsize=[6, 6])

    ax = fig.add_subplot(111)
    sigf = 2
    ax.set_title('$n_e = $' + tp(dens, sigf) + ' m$^{-3}$, $T = $' + tp(temp, sigf) + ' eV \n $B = $' + tp(bfield, sigf) +
                 ' T , $\\theta = $' + tp(viewangle * 180 / np.pi, sigf) + ' deg')

    line_models = ['rosato', 'stehle', 'stehle param', 'voigt']
    print('--------\ntimings:\n--------')

    for line_model in line_models:
        try:
            start_time = time.time()
            bls = pystark.BalmerLineshape(n_upper, dens, temp, bfield, viewangle=viewangle, line_model=line_model, wl_axis=wl_axis)
            end_time = time.time()

            print(line_model + ': ' + tp(end_time - start_time, sigf) + ' sec')
            ax.plot(wl_axis_nm, pystark.ls_norm(bls.ls_szd, wl_axis, norm_type=norm_type), '-', label=line_model)
        except Exception:
            print('--PYSTARK-- {} calculation failed:'.format(line_model))
            traceback.print_last()


    ax.set_xlim([np.min(wl_axis_nm), np.max(wl_axis_nm)])
    ax.axvline(pystark.get_wl_centre(n_upper) * 1e9, ls='--', color='dimgrey', zorder=0)
    leg = ax.legend(fontsize=fsize)
    ax.set_xlabel('wavelength (nm)', size=fsize)
    ax.set_yticklabels([])
    ax.set_yticks([])
    plt.semilogy()
    plt.show()

    return

if __name__ == '__main__':
    demo()
