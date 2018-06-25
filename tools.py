import numpy as np
import matplotlib.pyplot as plt

def find_nearest_idx(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


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

def norm(ls):
    return ls / np.max(ls)

