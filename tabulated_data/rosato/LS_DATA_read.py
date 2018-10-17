import numpy as np
import matplotlib.pyplot as plt
import pystark

import sys
import os
import time
sys.path.insert(0, pystark.paths.rosato_path)

import LS_DATA_read_copy

start = time.time()


# define inputs
n_upper = 5
dens = 1e14
temp = 1.25
bfield = 1.1
viewangle = 30

wmax = 5e-3
npts = 1001

detunings_axis = np.linspace(-wmax, wmax, npts)

# overwrite the in.txt fortran input file
LS_DATA_read_copy.in_param(n_upper, dens, temp, bfield, viewangle, wmax, npts)

# get the parameter grid bound indices
iN,iT,iB = LS_DATA_read_copy.set_bounds(dens, temp, bfield)

viewangle_idxs = [0, 1]
viewangle_rad = viewangle * np.pi / 180
ls = np.zeros([npts, 2])

for i, viewangle_idx in enumerate(viewangle_idxs):

    dir, name = LS_DATA_read_copy.set_name_file(n_upper, iN, iT, iB, viewangle_idx)

    w_arr, ls_arr = LS_DATA_read_copy.read_file(dir, name)
    ls[:, i] = LS_DATA_read_copy.ls_interpol(n_upper, dens, temp, bfield, wmax, npts, w_arr, ls_arr, iN, iT, iB)


ls = ls[:, 0] * np.sin(viewangle_rad) ** 2 + ls[:, 1] * np.cos(viewangle_rad) ** 2

end = time.time()

print(end - start, 'seconds')

plt.figure()
plt.plot(detunings_axis, ls)
plt.show()


