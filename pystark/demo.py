import pystark
from pystark.tools import to_precision as tp

import numpy as np
import os, sys, traceback, time

from scipy.constants import c
import matplotlib.pyplot as plt


def demo(wl_centre=None):
    """ Demo script compares all available lineshape models for selected Balmer line. """

    # specify plasma
    temp = 2  # [eV]
    dens = 1e20  # [m-3]
    bfield = 0.  # [T]
    viewangle = 90 * np.pi / 180  # [rad]

    # specify line
    isotope = 'D'
    n_upper = 5

    # plot params
    norm_type = 'area'
    fsize = 14

    # generate appropriate wavelength axis
    if wl_centre is None:
        wl_centre = pystark.get_wl_centre(n_upper)
    fwhm_voigt_hz = pystark.estimate_fwhm(n_upper, dens, temp, bfield, isotope)
    fwhm_voigt_m = c * fwhm_voigt_hz / (c / wl_centre) ** 2

    wl_axis = pystark.get_wl_axis(n_upper, dens, temp, bfield)
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
            bls = pystark.BalmerLineshape(n_upper, dens, temp, bfield, viewangle=viewangle, line_model=line_model, wl_axis=wl_axis, wl_centre=wl_centre)
            end_time = time.time()

            print(line_model + ': ' + tp(end_time - start_time, sigf) + ' sec')
            ax.plot(wl_axis_nm, pystark.ls_norm(bls.ls_szd, wl_axis, norm_type=norm_type), '-', label=line_model)
        except Exception:
            print('--PYSTARK-- {} calculation failed:'.format(line_model))
            traceback.print_last()

    ax.set_xlim([np.min(wl_axis_nm), np.max(wl_axis_nm)])
    ax.axvline(wl_centre * 1e9, ls='--', color='dimgrey', zorder=0)
    leg = ax.legend(fontsize=fsize)
    ax.set_xlabel('wavelength (nm)', size=fsize)
    ax.set_yticklabels([])
    ax.set_yticks([])
    # plt.semilogy()
    plt.show()

    return

if __name__ == '__main__':
    demo()
