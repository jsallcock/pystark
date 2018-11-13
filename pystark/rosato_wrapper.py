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


def rosato_wrapper(n_upper, dens, temp, bfield, viewangle, wmax, npts, display=False):
    """
     Stark-Zeeman lineshape profile for the first Balmer lines, interpolated from the Rosato et al. tabulated_data.
    
    Essentially a python wrapper for the LS_READ_data.f90
    
    :param n_upper: upper principal quantum number
    :param dens: [m^-3]
    :param temp: [eV]
    :param bfield: [T]
    :param viewangle: [rad] 
    :param wmax: [eV]
    :param npts: 
    :param display: bool.
    
    :return: 
    """

    dens_cm3 = dens / 1e6
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
    iN, iT, iB = LS_DATA_read_f2py.set_bounds(dens_cm3, temp, bfield)

    viewangle_idxs = [0, 1]
    ls = np.zeros([npts, 2])

    for i, viewangle_idx in enumerate(viewangle_idxs):
        name = LS_DATA_read_f2py.set_name_file(n_upper, iN, iT, iB, viewangle_idx)

        # replace the original 'dir' outputted from LS_DATA_read with one that is defined using the full path
        # -- this step could cause headaches later on, but allows fn to be called from anywhere.

        dir = os.path.join(rosato_database_path, balmer_line_names[n_upper] + '/').encode('utf-8')

        w_arr, ls_arr = LS_DATA_read_f2py.read_file(len(dir), dir, name)
        ls[:, i] = LS_DATA_read_f2py.ls_interpol(n_upper, dens_cm3, temp, bfield, wmax, npts, w_arr, ls_arr, iN, iT, iB)

    ls = ls[:, 0] * np.sin(viewangle) ** 2 + ls[:, 1] * np.cos(viewangle) ** 2

    if display:
        plt.figure()
        plt.plot(detunings_axis, ls)
        plt.ylabel('ls')
        plt.xlabel('Detuning (eV)')
        # plt.semilogy()
        plt.show()

    return detunings_axis, ls