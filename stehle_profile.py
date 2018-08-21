import numpy as np
import time
import pystark
from scipy.constants import *


def stehle_profile(n_upper, n_lower, temperature, density):
    """  Stark-broadened lineshape.
    
    :param n_upper: upper principal quantum number
    :param n_lower: lower principal quantum number
    :param temperature: in eV
    :param density: in m ** -3
    
    :return: wl_axis, lineshape_m, wprofs_nu2
    """

    # ensure given n_upper + n_lower fall within tabulated values
    assert n_lower in range(1, 4)
    assert n_upper in range(n_lower + 1, 31)

    temp_k = temperature * e / k  # temperature in K
    dens_cm = density * 1.e-6  # electronic density in cm-3
    prefix = 'n_' + str(n_upper) + '_' + str(n_lower) + '_'

    # extract raw tabulated data
    tab_temp_k = np.array(pystark.nc.variables[prefix + 'tempe'].data)  # tabulated electron temperatures (K)
    olam0 = pystark.nc.variables[prefix + 'olam0'].data  # line centre wavelength (A)
    num_tab_dens = pystark.nc.variables[prefix + 'id_max'].data
    fainom = pystark.nc.variables[prefix + 'fainom'].data
    tab_dens_cm = np.array(pystark.nc.variables[prefix + 'dense'].data)  # tabulated electron densities  (cm ** -3)
    f00 = np.array(pystark.nc.variables[prefix + 'f00'].data)      # normal Holtsmark field strength (30 kV / m)
    dl12 = np.array(pystark.nc.variables[prefix + 'dl12'].data)
    dl12s = np.array(pystark.nc.variables[prefix + 'dl12s'].data)
    fainu = pystark.nc.variables[prefix + 'fainu'].data  # Asymptotic value of iStark * (alpha ** 2.5) ("wings factor in alfa units")
    pr0 = np.array(pystark.nc.variables[prefix + 'pr0'].data)      # Ratio of the mean interelectronic distance to the electronic Debye length
    jtot = np.array(pystark.nc.variables[prefix + 'jtot'].data, dtype=np.int)  #  "number of wave lengths for the couple (T,Ne)"
    dom = np.array(pystark.nc.variables[prefix + 'dom'].data)  # frequency detunings in units (rad / (s*ues)
    d1om = np.array(pystark.nc.variables[prefix + 'd1om'].data)
    o1line = np.array(pystark.nc.variables[prefix + 'o1line'].data)
    o1lines = np.array(pystark.nc.variables[prefix + 'o1lines'].data)

    # ensure given temperature + density falls within tabulated values
    # change sligtly the value of the input density
    # dens_cm in order to remain , as far as possible, inside the tabulation
    # JSA: this first step seems bogus!

    if np.abs(dens_cm - tab_dens_cm[0]) / dens_cm <= 1.0E-3:
        dens_cm = tab_dens_cm[0] * 1.001

    for id in np.arange(1, num_tab_dens + 1):
        if np.abs(dens_cm - tab_dens_cm[id]) / dens_cm <= 1.0E-3:
            dens_cm = tab_dens_cm[id] * 0.999

    if dens_cm >= 2.0 * tab_dens_cm[num_tab_dens]:
        raise Exception('Your input density is higher than the largest tabulated value %f' % tab_dens_cm[num_tab_dens])

    if dens_cm <= tab_dens_cm[0]:
        raise Exception('Your input density is smaller than the smallest tabulated value %f' % tab_dens_cm[0])

    if temp_k >= tab_temp_k[9]:
        raise Exception('Your input temperature is higher than the largest tabulated value %f' % tab_temp_k[9])

    if temp_k <= tab_temp_k[0]:
        raise Exception('Your input temperature is lower than the smallest tabulated value %f' % tab_temp_k[0])

    normal_holtsmark_field = 1.25e-9 * (dens_cm ** (2. / 3.))  # normal field value in ues

    # calculate line centre wavelength and frequency using Rydberg formula
    # JSA: I have made this step clearer and corrected for deuteron mass in the Rydberg constant (though the effect is small)
    # TODO make sure this matches olam0 parameter above -- why were there two variables in the first place?!
    rydberg_m = Rydberg / (1. + (electron_mass / physical_constants['deuteron mass'][0]))
    wl_0_angst = 1e10 * (rydberg_m * (1 / n_lower ** 2 - 1 / n_upper ** 2)) ** -1

    # wl_0_angst = pystark.tools.get_NIST_balmer_wavelength(n_upper) * 1e10

    c_angst = c * 1e10  # velocity of light in Ansgtroms / s
    angular_freq_0 = 2 * np.pi * c_angst/ wl_0_angst  # rad / s

    otrans = -2 * np.pi * c_angst / wl_0_angst ** 2

    olines = o1lines / np.abs(otrans)
    oline = o1line / np.abs(otrans)

    # Limit analysis to uncorrelated plasmas.
    # check that mean interelectronic distance is smaller than the electronic Debye length (equ. 10)
    PR0_exp = 0.0898 * (dens_cm ** (1. / 6.)) / np.sqrt(temp_k)  # = (r0 / debye)
    if PR0_exp > 1.:
        raise Exception('The plasma is too strongly correlated\ni.e. r0/debye=0.1\nthe line cannot be computed.')

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

    # ==============================================
    # define an unique detunings grid - domm -  for the tabulated
    # profiles ( various temperatures , densities)
    # calculate all the line shapes for this  common grid
    # units used at this points are Domega_new= Delta(omega)/F00
    #                                      in rd/(s-1 ues)

    max_num_dens = 30  # Maximum number of densities
    max_num_tab_temp = 10
    max_num_detunings = 60  # Maximum number of detunings
    jtot = jtot.astype(np.int)
    domm = np.zeros(100000)
    dom0 = np.zeros(10000)
    tprof = np.zeros([max_num_dens, max_num_tab_temp, 10000])
    tprofs = np.zeros([max_num_dens, max_num_tab_temp, 10000])
    uprof = np.zeros([max_num_dens, 10000])
    uprofs = np.zeros([max_num_dens, 10000])

    inc = 0
    domm[inc] = 0.0
    # ---- Look to replace this loop
    for id in np.arange(num_tab_dens + 1):  # loop over tab densities
        for j in np.arange(max_num_tab_temp):  # loop over tab temperatures (?)
            for i in np.arange(1, jtot[id, j]):
                inc += 1
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
        dif = (dom0[i] - dom0[i - 1])
        if dif <= 1.0E-6:
            continue
        if dif / np.abs(dom0[i]) <= 0.1:
            continue
        inc = inc + 1
        domm[inc] = dom0[i]

    jdom = inc + 1  # One line after marker 35

    for id in np.arange(num_tab_dens):
        for j in np.arange(10):
            if pr0[id, j] > 1.0:
                continue

            tprof[id, j, 0] = oline[id, j, 0]
            tprofs[id, j, 0] = olines[id, j, 0]

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


    for id in np.arange(num_tab_dens):
        otest_dens = (dens_cm - tab_dens_cm[id]) * (dens_cm - tab_dens_cm[id + 1])
        if otest_dens <= 0.0:
            dense1 = tab_dens_cm[id]
            dense2 = tab_dens_cm[id + 1]
            id1 = id
            id2 = id + 1
            break

    if dens_cm >= tab_dens_cm[num_tab_dens]:
        dense1 = tab_dens_cm[num_tab_dens - 1]
        dense2 = tab_dens_cm[num_tab_dens]
        id1 = num_tab_dens - 1
        id2 = num_tab_dens

    for it in np.arange(10):
        otest = (temp_k - tab_temp_k[it]) * (temp_k - tab_temp_k[it + 1])
        if otest <= 0.0:
            it1 = it
            it2 = it + 1
            # pr01 = pr0[id2,it1] # max value of pr0 for T1,T2,dense1,dense2
            tempe1 = tab_temp_k[it]
            tempe2 = tab_temp_k[it + 1]
            break

    # interpolation in temperature
    for id in np.arange(id1, id2 + 1):
        for i in np.arange(jdom):
            uprof[id, i] = tprof[id, it1, i] + (temp_k - tempe1) * (tprof[id, it2, i] - tprof[id, it1, i]) / (
            tempe2 - tempe1)
            uprofs[id, i] = tprofs[id, it1, i] + (temp_k - tempe1) * (tprofs[id, it2, i] - tprofs[id, it1, i]) / (
            tempe2 - tempe1)

    delta_lambda = np.zeros(jdom)
    delta_nu = np.zeros(jdom)
    wprof_nu = np.zeros(jdom)
    wprofs_nu = np.zeros(jdom)

    for i in np.arange(jdom):
        wprof = uprof[id1, i] + (dens_cm - dense1) * (uprof[id2, i] - uprof[id1, i]) / (dense2 - dense1)
        wprofs = uprofs[id1, i] + (dens_cm - dense1) * (uprofs[id2, i] - uprofs[id1, i]) / (dense2 - dense1)
        delta_omega = domm[i] * normal_holtsmark_field
        delta_nu[i] = delta_omega / (2 * np.pi)
        delta_lambda[i] = wl_0_angst * delta_omega / (angular_freq_0 + delta_omega)
        # print(delta_lambda[i])
        wprof_nu[i] = (wprof / normal_holtsmark_field) * (2. * np.pi)
        wprofs_nu[i] = (wprofs / normal_holtsmark_field) * (2. * np.pi)
        #        print '%e %e %e %e' %(delta_lambda[i],delta_nu[i],wprof_nu[i],wprofs_nu[i])

    delta_lambda2 = np.concatenate((-delta_lambda[::-1], delta_lambda)) + olam0
    #    delta_nu2 = np.concatenate((-delta_nu[::-1],delta_nu))
    wprof_nu2 = np.concatenate((wprof_nu[::-1], wprof_nu))
    wprofs_nu2 = np.concatenate((wprofs_nu[::-1], wprofs_nu))

    # end = time.time()
    # print('Time elapsed (s): ', end - start)

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
