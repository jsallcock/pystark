import numpy as np
import matplotlib.pyplot as plt
import time
import pystark




def compare_ls(display=True):

    ls_path = '/Users/jsallcock/Documents/physics/phd/code/py_repo/pystark/tabulated_data/rosato/ls.txt'

    arr = np.loadtxt(ls_path)
    detunings_fortran = arr[:, 0]  # [eV]
    ls_fortran = arr[:, 1]  # looks like lineshape is area normalised to 1, although the wings are cut off for some grid-points.


    n_upper = 5
    dens = 2.41e13
    temp = 1.22
    bf = 0.026
    viewangle = 16.2
    wmax = 0.5e-3
    npts = 2000

    start = time.time()
    detunings_py, ls_py = pystark.rosato_stark_zeeman_profile(n_upper, dens, temp, bf, viewangle, wmax, npts, display=True)
    end = time.time()
    print(end - start, ' seconds')

    plt.figure()

    plt.plot(detunings_fortran, ls_fortran, '.-',  label='fortran')
    plt.plot(detunings_py, ls_py, '.-', label='python')
    plt.plot()
    plt.ylabel('ls')
    plt.xlabel('Detuning (eV)')


    plt.legend()
    plt.show()




    return


if __name__ == '__main__':
    compare_ls()