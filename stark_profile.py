import scipy.io.netcdf as netcdf
from scipy.constants import c, e, k
import matplotlib.pyplot as plt
import numpy as np
import time
import pystarky
from scipy.constants import *
import os
import scipy.integrate
from matplotlib.colors import LogNorm

# load data from 'stehle_tables.nc'
nc = netcdf.netcdf_file(pystarky.paths.netcdf_data_path, 'r')


def comparison():
    args = (4, 2, 1., 10.e20)

    x, y, ys = stehle_profile(*args)

    # get_fwhm(x, y, disp=True)

    plt.figure()

    plt.title('$n_e = $' + str(args[3]) + ' m$^{-3}$ \n $T = $' + str(args[2]) + ' eV')
    plt.plot(x, norm(y), '.-', label='Stehle')

    try:
        wl_axis, wl_0, lineshape_griem, lineshape_griem_stark, lineshape_griem_doppler = stark_profile(*args, mode='griem')

        plt.plot(wl_axis, norm(lineshape_griem), lw=1.25, label='Griem-Doppler profile (Voigt)')
        plt.plot(wl_axis, norm(lineshape_griem_stark), lw=0.75, label='Griem profile (Lorentzian)')
        plt.plot(wl_axis, norm(lineshape_griem_doppler), lw=0.75, label='Doppler profile (Gaussian)')
    except:
        print('Griem calculation failed')

    try:
        wl_axis, wl_0, lineshape_loman, lineshape_loman_stark, lineshape_loman_doppler = stark_profile(*args, mode='lomanowski')
        plt.plot(wl_axis, norm(lineshape_loman), lw=1.25, label='Lomanowski-Doppler profile (Voigt)')
        plt.plot(wl_axis, norm(lineshape_loman_stark), lw=0.75, label='Lomanowski profile (Lorentzian)')
    except:
        print('Lomanowski calculation failed')



    # plt.xlim(olam0 + np.array([-30.0, 30.0]))

    # plt.semilogy()
    # plt.xlim([np.min(wl_axis), np.max(wl_axis)])
    plt.legend()
    plt.xlabel('wavelength (m)')
    plt.show()


    return


def griem_error_load_data():

    import matplotlib.colors as colors

    loadpath = '/Users/jsallcock/Documents/physics/phd/code/py_repo/pystarky/saved_data/'

    err = np.load(os.path.join(loadpath, 'griem_error_percent.npy'))
    dens_axis = np.load(os.path.join(loadpath, 'griem_error_dens_axis_m3.npy'))
    temp_axis = np.load(os.path.join(loadpath, 'griem_error_temp_axis_ev.npy'))



    fig = plt.figure()
    ax = fig.add_subplot(111)
    # im = ax.pcolor(temp_axis, dens_axis, err, norm=LogNorm(vmin=np.nanmin(err), vmax=np.nanmax(err)))
    # cbar = fig.colorbar(im, ax=ax)

    cs = ax.contourf(temp_axis, dens_axis, err, [0.001, 0.01, 0.1, 1., 10, 100], norm=colors.LogNorm(vmin=np.nanmin(err), vmax=np.nanmax(err)))

    plt.colorbar(cs)
    # manual_locations = [(1., 1.)]
    # plt.clabel(cs, inline=1, fontsize=10)

    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.tight_layout()
    plt.show()


    return



# import matplotlib.cm as cm
#
# X = 10*np.random.rand(5,3)
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.imshow(X, cmap=cm.jet, interpolation='nearest')
#


# ax.format_coord = format_coord
# plt.show()


def griem_error_save_data():
    """ Compare the FWHM predicted by Griem's approximation to the Stehle tables. """

    # use full range of Stehle validity for comparison

    density_lims = (1e19, 0.9e25)  # m ** -3
    temperature_lims = (0.22, 105)  # eV

    n_upper, n_lower = (4, 2)

    len_dens_axis, len_temp_axis = (50, 40)
    dens_axis = np.logspace(np.log10(density_lims[0]), np.log10(density_lims[1]), len_dens_axis)
    temp_axis = np.logspace(np.log10(temperature_lims[0]), np.log10(temperature_lims[1]), len_temp_axis)

    err = np.zeros([len_dens_axis, len_temp_axis])

    for i, dens in enumerate(dens_axis):
        print(str(i) + ' / ' + str(len_dens_axis))
        for j, temp in enumerate(temp_axis):

            wl_axis, lineshape_griem, _, _ = stark_profile(n_upper, n_lower, temp, dens)
            try:
                x, y, ys = stehle_profile(n_upper, n_lower, temp, dens)
                fwhm_stehle = get_fwhm(x, y)
                fwhm_griem = get_fwhm(wl_axis, lineshape_griem)

                err[i, j] = 100 * abs(fwhm_stehle - fwhm_griem) / fwhm_stehle
            except:
                err[i, j] = np.nan


    savepath = '/Users/jsallcock/Documents/physics/phd/code/py_repo/pystarky/saved_data/'

    np.save(os.path.join(savepath, 'griem_error_percent.npy'), err)
    np.save(os.path.join(savepath, 'griem_error_dens_axis_m3.npy'), dens_axis)
    np.save(os.path.join(savepath, 'griem_error_temp_axis_ev.npy'), temp_axis)


    plt.figure()
    plt.imshow(err)
    plt.colorbar()
    plt.show()


    return


def find_nearest_idx(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


def get_fwhm(x, y, disp=False):
    """ given a function with a SINGLE PEAK, find the FWHM. QUICK AND DIRTY. """

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


def stark_profile(n_upper, n_lower, temperature, density, spectrum='wavelength', mode='griem'):
    """ Use Griem's scaling for an approximate Stark-broadened line in a certain parameter region. 
    
    - Scaling is well fulfilled for n_e > 10 ** 20 m ** -3 and T_i ~ 1eV.
    - Unverified for n_e > 10 ** 21 m ** -3
    
    TODO: CORRECT NORMALISATION TO MATCH STEHLE
    
    :param n_upper: upper principal quantum number
    :param n_lower: lower principal quantum number
    :param temperature: in eV
    :param density: in m ** -3
    
    """

    griem_valid_n_upper = [4, 5, 6, 8, 10]
    griem_valid_n_lower = [3]

    # lomanowski_valid_n_upper = []
    assert n_upper in griem_valid_n_upper

    valid_modes = ['griem', 'lomanowski']
    assert mode in valid_modes

    # load line centre from Stehle tables
    prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'
    wl_0 = nc.variables[prefix + 'olam0'].data[0] * 1e-10  # line centre wavlength (m)
    wl_0_nm = wl_0 * 1e9
    freq_0 = c / wl_0

    # ----- approximate FWHM of Voigt profile to dynamically generate frequency axis

    # Doppler FWHM: Gaussian
    temperature_k = temperature * e / k
    mass_kg = physical_constants['deuteron mass'][0]
    v_th = np.sqrt(2 * k * temperature_k / mass_kg)
    sigma_gauss_hz = v_th * freq_0 / (np.sqrt(2) * c)
    fwhm_gauss_hz = 2 * np.sqrt(2 * np.log(2)) * sigma_gauss_hz

    # Stark FWHM: Lorentzian
    griem_alpha_12s = {4: 0.08, 5: 0.09, 6: 0.17, 8: 0.28, 10: 0.46}
    alpha_12 = griem_alpha_12s[n_upper]

    ne_20 = density * 1e-20  # convert into units: 10 ** 20 m ** -3
    fwhm_lorentz_m = 1e-9 * 0.54 * alpha_12 * ne_20 ** (2. / 3.)  # Griem scaling
    fwhm_lorentz_hz = fwhm_lorentz_m * freq_0 ** 2 / c

    # total FWHM: Voigt
    fwhm_voigt_hz = 0.5346 * fwhm_lorentz_hz + np.sqrt(0.2166 * fwhm_lorentz_hz ** 2 + fwhm_gauss_hz ** 2)

    # generate frequency and wavelength axes
    num_fwhm = 20
    num_points = 20001  # odd number of points such that nu_0 lies exactly at the array centre.
    min_freq, max_freq = (freq_0 - num_fwhm * fwhm_voigt_hz, freq_0 + num_fwhm * fwhm_voigt_hz)
    freq_axis = np.linspace(min_freq, max_freq, num_points)

    min_wl, max_wl = (c / max_freq, c / min_freq)
    wl_axis = np.linspace(min_wl, max_wl, num_points)
    wl_axis_nm = wl_axis * 1e9

    # Doppler lineshape
    doppler_lineshape_hz = (freq_0 ** -1) * ((np.pi * (v_th ** 2 / (c ** 2))) ** -0.5) * np.exp(-(((freq_axis - freq_0) ** 2) / (2 * sigma_gauss_hz ** 2)))

    # Stark lineshape
    if mode == 'griem':

        hwhm_lorentz_hz = fwhm_lorentz_hz / 2
        stark_lineshape_hz = (1 / np.pi) * hwhm_lorentz_hz / ((freq_axis - freq_0) ** 2 + hwhm_lorentz_hz ** 2)

    elif mode == 'lomanowski':

        # Paramaterised MMM Stark profile coefficients from Bart's paper
        loman_abc_ij_idx = {'32': 0,
                            '42': 1,
                            '52': 2,
                            '62': 3,
                            '72': 4,
                            '82': 5,
                            '92': 6,
                            '43': 7,
                            '53': 8,
                            '63': 9,
                            '73': 10,
                            '83': 11,
                            '93': 12}
        loman_c_ij = [3.710e-18,
                      8.425e-18,
                      1.310e-15,
                      3.954e-16,
                      6.258e-16,
                      7.378e-16,
                      8.947e-16,
                      1.330e-16,
                      6.640e-16,
                      2.481e-15,
                      3.270e-15,
                      4.343e-15,
                      5.588e-15]
        loman_a_ij = [0.7665,
                      0.7803,
                      0.6796,
                      0.7149,
                      0.7120,
                      0.7159,
                      0.7177,
                      0.7449,
                      0.7356,
                      0.7118,
                      0.7137,
                      0.7133,
                      0.7165]
        loman_b_ij = [0.064,
                      0.050,
                      0.030,
                      0.028,
                      0.029,
                      0.032,
                      0.033,
                      0.045,
                      0.044,
                      0.016,
                      0.029,
                      0.032,
                      0.033]

        ij_idx = loman_abc_ij_idx[str(n_upper) + str(n_lower)]
        c_ij = loman_c_ij[ij_idx]
        a_ij = loman_a_ij[ij_idx]
        b_ij = loman_b_ij[ij_idx]

        delta_lambda_12ij = c_ij * (density ** a_ij) / (temperature ** b_ij)  # nm

        stark_lineshape_m = 1 / (abs(wl_axis_nm - wl_0_nm) ** (5. / 2.) + (delta_lambda_12ij / 2) ** (5. / 2.))
        stark_lineshape_m /= np.trapz(stark_lineshape_m, wl_axis)

        stark_lineshape_hz = np.interp(freq_axis, c / wl_axis[::-1], stark_lineshape_m[::-1] * (wl_axis ** 2 / c))

    # Voigt lineshape
    voigt_lineshape_hz = np.convolve(doppler_lineshape_hz, stark_lineshape_hz, 'same')
    voigt_lineshape_hz /= np.trapz(voigt_lineshape_hz, freq_axis)  # normalise


    # return spectra in wavelength or frequency:
    assert spectrum in ['wavelength', 'frequency']

    if spectrum == 'wavelength':

        doppler_lineshape_m = np.interp(wl_axis, c / freq_axis[::-1], doppler_lineshape_hz[::-1] * freq_axis ** 2 / c)
        stark_lineshape_m = np.interp(wl_axis, c / freq_axis[::-1], stark_lineshape_hz[::-1] * freq_axis ** 2 / c)
        voigt_lineshape_m = np.interp(wl_axis, c / freq_axis[::-1], voigt_lineshape_hz[::-1] * freq_axis ** 2 / c)

        return wl_axis, wl_0, voigt_lineshape_m, stark_lineshape_m, doppler_lineshape_m

    elif spectrum == 'frequency':
        return freq_axis, freq_0, voigt_lineshape_hz, stark_lineshape_hz, doppler_lineshape_hz


def stehle_profile(n_upper, n_lower, temperature, density, spectrum='wavelength'):
    """  Stark-broadened lineshape.
    
    :param n_upper: upper principal quantum number
    :param n_lower: lower principal quantum number
    :param temperature: in eV
    :param density: in m ** -3
    
    :return: 
    """

    # print('pystarky')

    start = time.time()

    assert n_lower in range(1, 4)
    assert n_upper in range(n_lower + 1, 31)

    prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'

    # extract raw_data
    tempe = np.array(nc.variables[prefix + 'tempe'].data)  # electron temperature (K)
    olam0 = nc.variables[prefix + 'olam0'].data  # line centre wavlength (A)
    id_max = nc.variables[prefix + 'id_max'].data
    fainom = nc.variables[prefix + 'fainom'].data
    dense = np.array(nc.variables[prefix + 'dense'].data)  # electron density (cm ** -3)
    f00 = np.array(nc.variables[prefix + 'f00'].data)      #  normal Holtsmark field strength (30 kV / m)
    dl12 = np.array(nc.variables[prefix + 'dl12'].data)
    dl12s = np.array(nc.variables[prefix + 'dl12s'].data)
    fainu = nc.variables[prefix + 'fainu'].data  # Asymptotic value of iStark * (alpha ** 2.5) ("wings factor in alfa units")
    pr0 = np.array(nc.variables[prefix + 'pr0'].data)      # Ratio of the mean interelectronic distance to the electronic Debye length
    jtot = np.array(nc.variables[prefix + 'jtot'].data, dtype=np.int)
    dom = np.array(nc.variables[prefix + 'dom'].data)
    d1om = np.array(nc.variables[prefix + 'd1om'].data)
    o1line = np.array(nc.variables[prefix + 'o1line'].data)
    o1lines = np.array(nc.variables[prefix + 'o1lines'].data)

    load_time = time.time() - start
    # print('load_time:', load_time)

    id_maxi = 30  # Maximum number of densities
    max_d = 60  # Maximum number of detunings

    TEMP = temperature * e / k  # temperature in K
    DENS = density * 1.e-6  # electronic density in cm-3

    jtot = jtot.astype(np.int)

    domm = np.zeros(100000)
    dom0 = np.zeros(10000)
    tprof = np.zeros([id_maxi, 10, 10000])
    tprofs = np.zeros([id_maxi, 10, 10000])
    uprof = np.zeros([id_maxi, 10000])
    uprofs = np.zeros([id_maxi, 10000])

    cspeed = c * 1e10  # velocity of light in Ansgtroms/s
    cspeed_pi = 2 * np.pi * cspeed

    PR0_exp = 0.0898 * (DENS ** (1. / 6.)) / np.sqrt(TEMP)  # =(r0/debye)
    F00_exp = 1.25E-9 * (DENS ** (2. / 3.))  # normal field value in ues

    ON = n_lower
    ONP = n_upper
    # DNU=1.0/ON**2 - 1.0/ONP**2
    ambda = 911.7633455 * (ON * ONP) ** 2 / ((ONP - ON) * (ONP + ON))

    omega = cspeed_pi / ambda
    otrans = -cspeed_pi / (ambda * ambda)

    olines = o1lines / np.abs(otrans)
    oline = o1line / np.abs(otrans)

    if PR0_exp > 1.:
        raise Exception('The plasma is too strongly correlated\ni.e. r0/debye=0.1\nthe line cannot be computed presently')

    # fainom_exp=fainom*(F00_exp**1.5)
    # fainum_exp=fainom_exp/( (OPI*2.)**1.5)

    # ========================
    # TABULATION Format CDS
    #   si on veut ecrire
    #  n -np lambda0 kalpha Ne E0 T R0/Debye Dalpha iDoppler iStark

    # IN_cds= N+0.01
    # INP_cds = NP+0.01

    # ***********************************************************
    # Don't edit the CDS format...
    # ***********************************************************

    # Skipped the code in the IF statement starting at line 470, since it
    # isn't used, if (.FALSE.) ...

    # ========================================================
    # change sligtly the value of the input density
    # DENS in order to remain , as far as possible, inside the tabulation

    if np.abs(DENS - dense[0]) / DENS <= 1.0E-3:
        DENS = dense[0] * 1.001

    for id in np.arange(1, id_max + 1):
        if np.abs(DENS - dense[id]) / DENS <= 1.0E-3:
            DENS = dense[id] * 0.999

    # ==============================================
    # define an unique detunings grid - domm -  for the tabulated
    # profiles ( various temperatures , densities)
    # calculate all the line shapes for this  common grid
    # units used at this points are Domega_new= Delta(omega)/F00
    #                                      in rd/(s-1 ues)

    inc = 0
    domm[inc] = 0.0
    # ---- Look to replace this loop
    for id in np.arange(id_max + 1):
        for j in np.arange(10):
            # print 'jtot[id,j]',jtot[id,j]
            for i in np.arange(1, jtot[id, j]):
                inc = inc + 1
                dom0[inc] = dom[id, j, i]

    inc = np.count_nonzero(dom)
    npik = inc + 1
    # nut=10000

    # Calling numpy sort instead of piksrt
    tmp = np.sort(dom0[0:npik])
    dom0[0:npik] = tmp[0:npik]
    # dom0 seems to agree with the FORTRAN version

    inc = 0
    domm[0] = 0.0
    # print 'npik',npik
    # ---- Look to replace this loop
    for i in np.arange(1, npik):
        # print 'i',i
        dif = (dom0[i] - dom0[i - 1])
        if dif <= 1.0E-6:
            continue
        if dif / np.abs(dom0[i]) <= 0.1:
            continue
        inc = inc + 1
        domm[inc] = dom0[i]

    jdom = inc + 1  # One line after marker 35

    for id in np.arange(id_max):
        for j in np.arange(10):
            if pr0[id, j] > 1.0:
                continue

            tprof[id, j, 0] = oline[id, j, 0]
            tprofs[id, j, 0] = olines[id, j, 0]

            # print 'id,j',id,j
            # print 'tprof,tprofs:',tprof[id,j,0],tprofs[id,j,0]

            if jtot[id, j] == 0:
                continue

            for i in np.arange(1, jdom + 1):
                skip1 = False
                skip2 = False
                # print 'i',i
                domeg = domm[i]
                ij_max = jtot[id, j]
                # print 'domeg,ij_max',domeg,ij_max
                for ij in np.arange(1, ij_max - 1):
                    # print 'ij',ij
                    test = (domeg - dom[id, j, ij]) * (domeg - dom[id, j, ij - 1])
                    # print 'test1:',test
                    if test <= 0.0:
                        # print 'triggered test1'
                        x1 = dom[id, j, ij - 1]
                        x2 = dom[id, j, ij]
                        x3 = dom[id, j, ij + 1]
                        y1 = oline[id, j, ij - 1]
                        y2 = oline[id, j, ij]
                        y3 = oline[id, j, ij + 1]
                        # print 'x1,x2,x3',x1,x2,x3
                        # print 'y1,y2,y3',y1,y2,y3
                        tprof[id, j, i] = FINTRP(x1, x2, x3, y1, y2, y3, domeg)
                        y1 = olines[id, j, ij - 1]
                        y2 = olines[id, j, ij]
                        y3 = olines[id, j, ij + 1]
                        tprofs[id, j, i] = FINTRP(x1, x2, x3, y1, y2, y3, domeg)
                        # print 'tprof[id,j,i]',tprof[id,j,i]
                        # print 'tprofs[id,j,i]',tprofs[id,j,i]
                        skip1 = True
                        skip2 = True
                        break

                if skip1 is False:
                    test = (domeg - dom[id, j, ij_max - 2]) * (domeg - dom[id, j, ij_max - 1])
                    # print 'test2:',test
                    # print 'domeg',domeg
                    # print 'dom[id,j,ij_max-1]',dom[id,j,ij_max-2]
                    # print 'dom[id,j,ij_max]',dom[id,j,ij_max-1]
                    if test <= 0.0:
                        # print 'triggered test2'
                        x1 = dom[id, j, ij_max - 3]
                        x2 = dom[id, j, ij_max - 2]
                        x3 = dom[id, j, ij_max - 1]
                        y1 = oline[id, j, ij_max - 3]
                        y2 = oline[id, j, ij_max - 2]
                        y3 = oline[id, j, ij_max - 1]
                        tprof[id, j, i] = FINTRP(x1, x2, x3, y1, y2, y3, domeg)
                        y1 = olines[id, j, ij_max - 3]
                        y2 = olines[id, j, ij_max - 2]
                        y3 = olines[id, j, ij_max - 1]
                        tprofs[id, j, i] = FINTRP(x1, x2, x3, y1, y2, y3, domeg)
                        skip2 = True
                        # print 'x1,x2,x3',x1,x2,x3
                        # print 'y1,y2,y3',y1,y2,y3
                        # print 'tprof[id,j,i]',tprof[id,j,i]
                        # print 'tprofs[id,j,i]',tprofs[id,j,i]
                        continue

                if skip2 is False:
                    if domeg > dom[id, j, ij_max]:
                        # print 'triggered test3'
                        tprof[id, j, i] = fainom / (domeg ** 2.5)
                        tprofs[id, j, i] = tprof[id, j, i]
                        continue

    # We can skip writing the intermediate file

    if DENS >= 2.0 * dense[id_max]:
        raise Exception('Your input density is higher than the largest tabulated value %f' % dense[id_max])

    if DENS <= dense[0]:
        raise Exception('Your input density is smaller than the smallest tabulated value %f' % dense[0])

    if TEMP >= tempe[9]:
        raise Exception('Your input temperature is higher than the largest tabulated value %f' % tempe[9])

    if TEMP <= tempe[0]:
        raise Exception('Your input temperature is lower than the smallest tabulated value %f' % tempe[0])

    for id in np.arange(id_max):
        otest_dens = (DENS - dense[id]) * (DENS - dense[id + 1])
        if otest_dens <= 0.0:
            dense1 = dense[id]
            dense2 = dense[id + 1]
            id1 = id
            id2 = id + 1
            break

    if DENS >= dense[id_max]:
        dense1 = dense[id_max - 1]
        dense2 = dense[id_max]
        id1 = id_max - 1
        id2 = id_max

    for it in np.arange(10):
        otest = (TEMP - tempe[it]) * (TEMP - tempe[it + 1])
        if otest <= 0.0:
            it1 = it
            it2 = it + 1
            # pr01 = pr0[id2,it1] # max value of pr0 for T1,T2,dense1,dense2
            tempe1 = tempe[it]
            tempe2 = tempe[it + 1]
            break

    # interpolation in temperature
    for id in np.arange(id1, id2 + 1):
        for i in np.arange(jdom):
            uprof[id, i] = tprof[id, it1, i] + (TEMP - tempe1) * (tprof[id, it2, i] - tprof[id, it1, i]) / (
            tempe2 - tempe1)
            uprofs[id, i] = tprofs[id, it1, i] + (TEMP - tempe1) * (tprofs[id, it2, i] - tprofs[id, it1, i]) / (
            tempe2 - tempe1)

    delta_lambda = np.zeros(jdom)
    delta_nu = np.zeros(jdom)
    wprof_nu = np.zeros(jdom)
    wprofs_nu = np.zeros(jdom)

    for i in np.arange(jdom):
        wprof = uprof[id1, i] + (DENS - dense1) * (uprof[id2, i] - uprof[id1, i]) / (dense2 - dense1)
        wprofs = uprofs[id1, i] + (DENS - dense1) * (uprofs[id2, i] - uprofs[id1, i]) / (dense2 - dense1)
        delta_omega = domm[i] * F00_exp
        delta_nu[i] = delta_omega / (2 * np.pi)
        delta_lambda[i] = ambda * delta_omega / (omega + delta_omega)
        # print(delta_lambda[i])
        wprof_nu[i] = (wprof / F00_exp) * (2. * np.pi)
        wprofs_nu[i] = (wprofs / F00_exp) * (2. * np.pi)
        #        print '%e %e %e %e' %(delta_lambda[i],delta_nu[i],wprof_nu[i],wprofs_nu[i])

    delta_lambda2 = np.concatenate((-delta_lambda[::-1], delta_lambda)) + olam0
    #    delta_nu2 = np.concatenate((-delta_nu[::-1],delta_nu))
    wprof_nu2 = np.concatenate((wprof_nu[::-1], wprof_nu))
    wprofs_nu2 = np.concatenate((wprofs_nu[::-1], wprofs_nu))

    end = time.time()
    print('Time elapsed (s): ', end - start)

    # area normalise and convert wavelength units to m:
    wl_axis = delta_lambda2 * 1e-10
    lineshape_m = wprof_nu2 / np.trapz(wprof_nu2, wl_axis)


    return wl_axis, lineshape_m, wprofs_nu2


def FINTRP(x1, x2, x3, y1, y2, y3, x):

    if x == x2:
        return y2

    a12 = x1 - x2
    a22 = x1 - x3
    v1 = y1 - y2
    v2 = y1 - y3

    if ((y1 < y2) and (y2 < y3)) or ((y1 > y2) and (y2 > y3)):
        deter = v1 * a22 - v2 * a12
        if np.abs(deter) < 1.0E-40:
            return y1 + (x - x1) * (y3 - y1) / (x3 - x1)
        a21 = x1 * y1
        a11 = a21 - x2 * y2
        a21 = a21 - x3 * y3
        c = (a22 * a11 - a12 * a21) / deter
        a = (-v2 * a11 + v1 * a21) / deter
        b = (y1 - a) * (x1 - c)
        return a + b / (x - c)
    else:
        x1c = x1 * x1
        a11 = x1c - x2 * x2
        a21 = x1c - x3 * x3
        deter = a11 * a22 - a12 * a21
        if np.abs(deter) < 1.0E-40:
            raise Exception('FINTRP error: incorrect inputs')
        a = (a22 * v1 - a12 * v2) / deter
        b = (-a21 * v1 + a11 * v2) / deter
        c = y1 - a * x1c - b * x1
        return (a * x + b) * x + c


def norm(ls):
    return ls / np.max(ls)


if __name__ == '__main__':
    comparison()
    # griem_error_save_data()
    # griem_error_load_data()