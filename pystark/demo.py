import matplotlib
gui_env = ['MacOSX', 'TKAgg', 'GTKAgg', 'Qt4Agg', 'WXAgg']
for gui in gui_env:
    try:
        matplotlib.use(gui, warn=False, force=True)
        from matplotlib import pyplot as plt
        break
    except:
        continue
from pystark.tools import to_precision as tp
import pystark
import numpy as np
import os, sys, traceback, time

from scipy.constants import c
import matplotlib.pyplot as plt


def demo(wl_centre=None):
    """ Demo script compares all available lineshape models for selected Balmer line. """

    # specify plasma
    temp = 2.  # [eV]
    dens = 8e20  # [m-3]
    bfield = 0.  # [T]
    viewangle = 90 * np.pi / 180  # [rad]

    # specify line
    species = 'H'
    n_upper = 5
    n_lower = 2

    # plot params
    norm_type = 'area'
    fsize = 10

    # generate appropriate wavelength axis
    if wl_centre is None:
        wl_centre = pystark.get_wl_centre(species, n_upper, n_lower)
    wl_axis = pystark.get_wl_axis(species, n_upper, n_lower, dens, temp, bfield)
    wl_axis_nm = wl_axis * 1e9

    # generate plots
    fig = plt.figure(figsize=[10, 5])

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    sigf = 2
    fig.suptitle('$n_e = $' + tp(dens, sigf) + ' m$^{-3}$, $T = $' + tp(temp, sigf) + ' eV \n $B = $' + tp(bfield, sigf) +
                 ' T , $\\theta = $' + tp(viewangle * 180 / np.pi, sigf) + ' deg')

    line_models = ['rosato', 'stehle', 'stehle_param', 'voigt']
    print('--------\ntimings:\n--------')

    for line_model in line_models:
        try:
            start_time = time.time()
            bls = pystark.StarkLineshape(species, n_upper, n_lower, dens, temp, bfield, view_angle=viewangle, line_model=line_model,
                                         wl_axis=wl_axis, wl_centre=wl_centre, override_input_check=True)
            end_time = time.time()

            print(line_model + ': ' + tp(end_time - start_time, sigf) + ' sec')
            ax1.plot(wl_axis_nm, pystark.ls_norm(bls.ls_szd, wl_axis, norm_type=norm_type), '-', label=line_model)
            ax2.plot(wl_axis_nm, pystark.ls_norm(bls.ls_szd, wl_axis, norm_type=norm_type), '-', label=line_model)
        except:
            print('--PYSTARK-- {} calculation failed:'.format(line_model))
            # traceback.print_last()
    axes = [ax1, ax2]
    for ax in axes:
        ax.set_xlim([np.min(wl_axis_nm), np.max(wl_axis_nm)])
        ax.axvline(wl_centre * 1e9, ls='--', color='dimgrey', zorder=0)
        leg = ax.legend(fontsize=fsize)
        ax.set_xlabel('wavelength (nm)', size=fsize)
        ax.set_yticklabels([])
        ax.set_yticks([])

    ax2.semilogy()
    # plt.semilogy()
    plt.show()

    return


def demo_helium():
    pass

    return


if __name__ == '__main__':
    demo()
