import pystark as ps
import numpy as np
import matplotlib.pyplot as plt


def plot_comparison():
    """ Demo script compares all lineshape models for selected Balmer line. """

    n_upper = 5
    n_lower = 2
    temp = 10.  # [eV]
    dens = 1e19  # [m-3]
    bfield = 0.  # [T]
    viewangle = 0  # [deg]

    # generate appropriate wavelength axis

    wl_0 = ps.get_NIST_balmer_wavelength(n_upper)
    fwhm_voigt_hz, _, _ = ps.estimate_fwhms(n_upper, dens, temp)
    fwhm_voigt_m = 3e8 * fwhm_voigt_hz / (3e8 / wl_0) ** 2

    wl_axis = np.linspace(wl_0 - 5 * fwhm_voigt_m, wl_0 + 5 * fwhm_voigt_m, 200)
    wl_axis_nm = wl_axis * 1e9

    display=False

    fig = plt.figure(figsize=[6, 6])
    fsize=14
    ax = fig.add_subplot(111)
    ax.set_title('$n_e = $' + str(dens) + ' m$^{-3}$ \n $T = $' + str(temp) + ' eV \n $B = $' + str(bfield) + ' T')

    try:
        ls_stehle = ps.stehle_profile(n_upper, n_lower, temp, dens, wl_axis, display=display)
        ax.plot(wl_axis_nm, ps.norm_a(ls_stehle, wl_axis), '-', label='Stehle')
    except:
        print('Stehle calculation failed')

    try:
        ls_rosato = ps.rosato_profile(n_upper, dens, temp, bfield, viewangle, wl_axis, display=display)
        ax.plot(wl_axis_nm, ps.norm_a(ls_rosato, wl_axis), '-', label='Rosato')
    except:
        print('Rosato calculation failed')

    try:
        ls_loman, _, _ = ps.simple_profile(n_upper, n_lower, temp, dens, wl_axis, model='lomanowski')
        ax.plot(wl_axis_nm, ps.norm_a(ls_loman, wl_axis), '-', label='Lomanowski')
    except:
        print('Lomanowski calculation failed')

    try:
        ls_griem, _, _ = ps.simple_profile(n_upper, n_lower, temp, dens, wl_axis, model='griem')
        ax.plot(wl_axis_nm, ps.norm_a(ls_griem, wl_axis), '-', label='Griem')
    except:
        print('Griem calculation failed')


    ax.set_xlim([np.min(wl_axis_nm), np.max(wl_axis_nm)])
    leg = ax.legend(fontsize=fsize)
    ax.set_xlabel('wavelength (nm)', size=fsize)
    # plt.semilogy()

    return


if __name__ == '__main__':
    plot_comparison()
    plt.show()
