import numpy as np
import matplotlib.pyplot as plt
import os
import pystark

# parameter grid, as defined at the start of section 3
dens_grid = np.array([1.e13, 2.15e13, 4.64e13, 1.e14, 2.15e14, 4.64e14, 1.e15, 2.15e15, 4.64e15, 1.e16])  # [cm-3]
temp_grid = np.array([0.316, 1., 3.16, 10, 31.6])  # [eV]
bfield_grid = np.array([0., 1., 2., 2.5, 3., 5.])  # [T]
viewangle_grid = np.array([90., 0.])  # [rad.]


def rosato_stark_zeeman_profile_old(n_upper, dens, temp, bf, viewangle, wmax, npts, display=False):
    """ NOW OBSOLETE, REPLACED BY rosato_wrapper.py, WHICH USE FORTRAN WRAPPER TO ACHIEVE >100X SPEEDUP 
    
    A python port of the LS_DATA_read.f90 from the Rosato et al. supplementary material. 

    Oh boy, this is mostly translated line-by-line, can you tell?  

    BENCHMARKED AGAINST THE FORTRAN CODE, ls and detunings values match fortran output to 1 part in ~10^6. Useful to 
    keep for comparison/ future benchmarking.

    :param n_upper: 
    :param dens: [cm-3]
    :param temp: [eV]
    :param bf: [T]
    :param viewangle: [deg]
    :param wmax: max detuning [eV]
    :param npts: length of detuning axis

    :return: 
    """

    assert dens > dens_grid[0] and dens < dens_grid[-1]
    assert temp > temp_grid[0] and temp < temp_grid[-1]
    assert bf >= bfield_grid[0] and bf < bfield_grid[-1]

    # user-defined detunings grid
    detuning_axis = np.linspace(-wmax, wmax, npts)

    w_arr = np.zeros([2, 2, 2, 2, 1000])

    ls_arr = np.zeros([2, 2, 2, 2, 1000])
    ls_arr2 = np.zeros([2, 2, 2, 2, npts])
    ls_arr3 = np.zeros([2, 2, 2, npts])
    ls_arr4 = np.zeros([2, npts])

    # find upper and lower bounding indices for the relevant parameter grids
    temp_1_idx = np.argmax(temp_grid > temp)
    temp_0_idx = temp_1_idx - 1
    temp_idxs = [temp_0_idx, temp_1_idx]

    dens_1_idx = np.argmax(dens_grid > dens)
    dens_0_idx = dens_1_idx - 1
    dens_idxs = [dens_0_idx, dens_1_idx]

    bfield_ub_idx = np.argmax(bfield_grid > bf)
    bfield_lb_idx = bfield_ub_idx - 1
    bfield_idxs = [bfield_lb_idx, bfield_ub_idx]

    # load all tabulated_data from .txt files necessary for the interpolation
    for i_angle in range(2):
        for i_dens in range(2):
            for i_temp in range(2):
                if i_dens * i_temp == 1:
                    break
                    # break ensures that only the necessary lineshape tabulated_data is loaded (i_dens = 1, i_temp = 1 is not used)
                for i_bfield in range(2):
                    w_arr[i_dens, i_temp, i_bfield, i_angle, :], ls_arr[i_dens, i_temp, i_bfield, i_angle, :] \
                        = read_single_file(n_upper, dens_idxs[i_dens], temp_idxs[i_temp], bfield_idxs[i_bfield], i_angle)



    # interpolation
    if bfield_ub_idx != 1:
        for i_angle in range(2):
            for i, w in enumerate(detuning_axis):
                for i_dens in range(2):
                    for i_temp in range(2):
                        if i_dens * i_temp == 1:
                            break
                            # break ensures that only the necessary lineshape tabulated_data is loaded (i_dens = 1, i_temp = 1 is not used)
                        for i_bfield in range(2):
                            for iw in range(1000):
                                if w < (w_arr[i_dens, i_temp, i_bfield, i_angle, iw] * bf / bfield_grid
                                            [bfield_ub_idx - 1 + i_bfield]):
                                    break  # ensures that the lineshape is zero outside of the tabulated detuning values

                            if iw != 0 and iw < 999:
                                ls_arr2[i_dens, i_temp, i_bfield, i_angle, i] = ls_arr[i_dens, i_temp, i_bfield, i_angle, iw -1] + \
                                                                                (w * bfield_grid[bfield_ub_idx - 1 + i_bfield] / bf - w_arr[i_dens, i_temp, i_bfield, i_angle, iw - 1]) * \
                                                                                (ls_arr[i_dens, i_temp, i_bfield, i_angle, iw] - ls_arr[i_dens, i_temp, i_bfield, i_angle, iw - 1]) / \
                                                                                (w_arr[i_dens, i_temp, i_bfield, i_angle, iw] - w_arr[i_dens, i_temp, i_bfield, i_angle, iw - 1])

                        # interpolation in B-field (equ. 3) TODO: CLEAN UP
                        ls_arr3[i_dens, i_temp, i_angle, i] = ((bf - bfield_grid[bfield_ub_idx - 1]) / (
                        bfield_grid[bfield_ub_idx] - bfield_grid[bfield_ub_idx - 1])) * \
                                                              (bfield_grid[bfield_ub_idx] / bf) * ls_arr2[
                                                                  i_dens, i_temp, 1, i_angle, i] + \
                                                              ((bfield_grid[bfield_ub_idx] - bf) / (
                                                              bfield_grid[bfield_ub_idx] - bfield_grid[
                                                                  bfield_ub_idx - 1])) * \
                                                              (bfield_grid[bfield_ub_idx - 1] / bf) * ls_arr2[
                                                                  i_dens, i_temp, 0, i_angle, i]

    elif bf == 0.:
        for i, w in enumerate(detuning_axis):
            for i_angle in range(2):
                for i_dens in range(2):
                    for i_temp in range(2):
                        if i_dens * i_temp == 1:
                            break
                            # break ensures that only the necessary lineshape tabulated_data is loaded (i_dens = 1, i_temp = 1 is not used)
                        for iw in range(1000):
                            if w < w_arr[i_dens, i_temp, 0, i_angle, iw]:
                                break  # ensures that the lineshape is zero outside of the tabulated detuning values

                        if iw != 0 and iw < 999:
                            ls_arr2[i_dens, i_temp, 0, i_angle, i] = ls_arr[i_dens, i_temp, 0, i_angle, iw - 1] + \
                                                                     (w - w_arr[i_dens, i_temp, 0, i_angle, iw - 1]) * \
                                                                     (ls_arr[i_dens, i_temp, 0, i_angle, iw] - ls_arr[
                                                                         i_dens, i_temp, 0, i_angle, iw - 1]) / \
                                                                     (w_arr[i_dens, i_temp, 0, i_angle, iw] - w_arr[
                                                                         i_dens, i_temp, 0, i_angle, iw - 1])

                        ls_arr3[i_dens, i_temp, i_angle, i] = ls_arr2[i_dens, i_temp, 0, i_angle, i]

    else:
        for i, w in enumerate(detuning_axis):
            for i_angle in range(2):
                for i_dens in range(2):
                    for i_temp in range(2):
                        if i_dens * i_temp == 1:
                            break
                            # break ensures that only the necessary lineshape tabulated_data is loaded (i_dens = 1, i_temp = 1 is not used)
                        for iw in range(1000):
                            if w < w_arr[i_dens, i_temp, 0, i_angle, iw]:
                                break  # ensures that the lineshape is zero outside of the tabulated detuning values

                        if iw != 0 and iw < 999:
                            ls_arr2[i_dens, i_temp, 0, i_angle, i] = ls_arr[i_dens, i_temp, 0, i_angle, iw - 1] + (w -
                                                                                                                   w_arr[
                                                                                                                       i_dens, i_temp, 0, i_angle, iw - 1]) * \
                                                                                                                  (
                                                                                                                  ls_arr[
                                                                                                                      i_dens, i_temp, 0, i_angle, iw] -
                                                                                                                  ls_arr[
                                                                                                                      i_dens, i_temp, 0, i_angle, iw - 1]) / \
                                                                                                                  (
                                                                                                                  w_arr[
                                                                                                                      i_dens, i_temp, 0, i_angle, iw] -
                                                                                                                  w_arr[
                                                                                                                      i_dens, i_temp, 0, i_angle, iw - 1])
                        for iw in range(1000):
                            if w < (w_arr[i_dens, i_temp, 1, i_angle, iw] * bf / bfield_grid[1]):
                                break

                        if iw != 0 and iw < 999:
                            ls_arr2[i_dens, i_temp, 1, i_angle, i] = ls_arr[i_dens, i_temp, 1, i_angle, iw - 1] + \
                                                                     (w * bfield_grid[1] / bf - w_arr[
                                                                         i_dens, i_temp, 1, i_angle, iw - 1]) * \
                                                                     (ls_arr[i_dens, i_temp, 1, i_angle, iw] - ls_arr[
                                                                         i_dens, i_temp, 1, i_angle, iw - 1]) / \
                                                                     (w_arr[i_dens, i_temp, 1, i_angle, iw] - w_arr[
                                                                         i_dens, i_temp, 1, i_angle, iw - 1])

                        ls_arr3[i_dens, i_temp, i_angle, i] = ((bf - bfield_grid[bfield_ub_idx - 1]) / (
                        bfield_grid[bfield_ub_idx] - bfield_grid[bfield_ub_idx - 1])) * \
                                                              (bfield_grid[bfield_ub_idx] / bf) * ls_arr2[
                                                                  i_dens, i_temp, 1, i_angle, i] + \
                                                              ((bfield_grid[bfield_ub_idx] - bf) / (
                                                              bfield_grid[bfield_ub_idx] - bfield_grid[
                                                                  bfield_ub_idx - 1])) * \
                                                              ls_arr2[i_dens, i_temp, 0, i_angle, i]

    # interpolation in density and temperature (equ. 2)
    ls_arr4[:, :] = 3. * np.log10(dens / dens_grid[dens_1_idx - 1]) * ls_arr3[1, 0, :, :] + \
                    2. * np.log10(temp / temp_grid[temp_1_idx - 1]) * ls_arr3[0, 1, :, :] + \
                    (1. - 3. * np.log10(dens / dens_grid[dens_1_idx - 1]) - 2. * np.log10(
                        temp / temp_grid[temp_1_idx - 1])) * ls_arr3[0, 0, :, :]

    # account for view angle (section 3)
    viewangle *= 3.141593 / 180.  # np.pi / 180  # [deg] to [rad]
    ls = ls_arr4[0, :] * np.sin(viewangle) ** 2 + ls_arr4[1, :] * np.cos(viewangle) ** 2

    if display:
        plt.figure()
        plt.plot(detuning_axis, ls)
        plt.ylabel('ls')
        plt.xlabel('Detuning (eV)')
        plt.show()

    return detuning_axis, ls


def read_single_file(n_upper, dens_idx, temp_idx, bfield_idx, viewangle_idx):
    """ Given the indices of the plasma parameter grids, load the lineshape from file.

    indexing inputs begins at 0 here, whereas in the filenames it begins at 1 ...except for viewangle...

    :param n_upper: integer within [3, 7]
    :param dens_idx: integer within [0, 9]
    :param temp_idx: integer within [0, 4]
    :param bfield_idx: integer within [0, 5]
    :param viewangle_idx: [0, 1]
    :return: 
    """

    balmer_line_names = ['-', '-', '-', 'D_alpha', 'D_beta', 'D_gamma', 'D_delta',
                         'D_epsilon']  # list index gives transition n_upper

    load_dir_path = os.path.join(pystark.paths.rosato_database_path, balmer_line_names[n_upper])
    filename = 'ls%02d' % (dens_idx + 1,) + str(temp_idx + 1) + str(bfield_idx + 1) + str(viewangle_idx) + '.txt'
    load_file_path = os.path.join(load_dir_path, filename)

    arr = np.loadtxt(load_file_path)

    return arr[:, 0], arr[:, 1]