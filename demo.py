import pystark as ps
import numpy as np
import matplotlib.pyplot as plt


def plot_comparison():
    """ Demo script, compare the Stehle lineshape with a Voigt profile (based on Griem's FWHM density scaling) and  """
    args = (4, 2, 1., 1.e20)  # (n_upper, n_lower, temperature [eV], Density [m ** -3])

    x, y, ys = ps.stehle_profile(*args)


    plt.figure()
    plt.title('$n_e = $' + str(args[3]) + ' m$^{-3}$ \n $T = $' + str(args[2]) + ' eV')
    plt.plot(x, ps.norm_h(y), '.-', label='Stehle')

    try:
        wl_axis, wl_0, lineshape_griem, lineshape_griem_stark, lineshape_griem_doppler = ps.simple_profile(*args, model='griem')

        plt.plot(wl_axis, ps.norm_h(lineshape_griem), lw=1.25, label='Griem-Doppler profile (Voigt)')
        plt.plot(wl_axis, ps.norm_h(lineshape_griem_stark), ':', lw=0.75, label='Griem profile (Lorentzian)')
        plt.plot(wl_axis, ps.norm_h(lineshape_griem_doppler), ':', lw=0.75, label='Doppler profile (Gaussian)')
    except:
        print('Griem calculation failed')

    try:
        wl_axis, wl_0, lineshape_loman, lineshape_loman_stark, lineshape_loman_doppler = ps.simple_profile(*args, model='lomanowski')
        plt.plot(wl_axis, ps.norm_h(lineshape_loman), lw=1.25, label='Lomanowski-Doppler profile (Voigt)')
        plt.plot(wl_axis, ps.norm_h(lineshape_loman_stark), ':', lw=0.75, label='Lomanowski profile (Lorentzian)')
    except:
        print('Lomanowski calculation failed')

    plt.xlim([np.min(wl_axis), np.max(wl_axis)])
    plt.legend()
    plt.xlabel('wavelength (m)')

    return


if __name__ == '__main__':
    plot_comparison()
    plt.show()
