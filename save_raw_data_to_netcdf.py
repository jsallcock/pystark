# -*- coding: utf-8 -*-\
"""
Created on Sat Feb 03 17:18:06 2018

@author: jrh

Edited by jsallcock 11/06/18
"""

import pystarky
import os, fnmatch
import re
import numpy as np
from scipy.io import netcdf

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def fort_conv(string):
    if bool(string.strip()):
        return np.float(string)
    else:
        return 0.0



write_netcdf = True

if write_netcdf:
    ncf = netcdf.netcdf_file(pystarky.paths.netcdf_data_path, 'w')
    ncf.createDimension('num', 1)
    ncf.createDimension('ten', 10)
    ncf.createDimension('max_d', 60)


# path details for finding directories
rdp = pystarky.paths.raw_data_path

ly_prefix, ly_n_l_lim, ly_n_u_lim = 'ly', 1, 30  # lyman
ba_prefix, ba_n_l_lim, ba_n_u_lim = 'ba', 2, 30  # balmer
pa_prefix, pa_n_l_lim, pa_n_u_lim = 'pa', 3, 30  # paschen

ly_num = ly_n_u_lim - ly_n_l_lim
ba_num = ba_n_u_lim - ba_n_l_lim
pa_num = pa_n_u_lim - pa_n_l_lim

transition_num = ly_num + ba_num + pa_num

for transition in range(0, transition_num):

    # specify directory
    if transition in range(0, ly_num):
        dir = os.path.join(rdp, ly_prefix, ly_prefix + "%02d" % (transition + 1 + ly_n_l_lim,))

    elif transition in range(ly_num, ly_num + ba_num):
        dir = os.path.join(rdp, ba_prefix, ba_prefix + "%02d" % (transition + 1 + ba_n_l_lim - ly_num,))

    elif transition in range(ly_num + ba_num, transition_num):
        dir = os.path.join(rdp, pa_prefix, pa_prefix + "%02d" % (transition + 1 + pa_n_l_lim - ly_num - ba_num,))

    id_maxi = len(find('profil*.dat', dir))  # Maximum number of densities
    max_d = 60  # Maximum number of detunings

    TEMP = 10000.0  # temperature in K
    DENS = 5.0E14  # electronic density in cm-3

    tempe = np.zeros(10)
    jtot = np.zeros([id_maxi,10])
    din = np.zeros([10,max_d])
    sprof = np.zeros([10,max_d])
    sprofs = np.zeros([10,max_d])
    dl12 = np.zeros(10)
    dl12s = np.zeros(10)

    dom = np.zeros([id_maxi,10,max_d])
    oline = np.zeros([id_maxi,10,max_d])
    olines = np.zeros([id_maxi,10,max_d])

    d1om = np.zeros([id_maxi,10,max_d])
    o1line = np.zeros([id_maxi,10,max_d])
    o1lines = np.zeros([id_maxi,10,max_d])

    dom0 = np.zeros(10000)
    domm = np.zeros(100000)
    tprof = np.zeros([id_maxi,10,10000])
    tprofs = np.zeros([id_maxi,10,10000])
    dense = np.zeros(id_maxi)
    f00 = np.zeros(id_maxi)
    pr0 = np.zeros([id_maxi,10])
    uprof = np.zeros([id_maxi,10000])
    uprofs = np.zeros([id_maxi,10000])

    with open(os.path.join(dir, 'nraie.dat')) as f:
        id_max = np.int(f.readline()) - 1  # number of density-points for the transition
        f.close()

    # OPI=3.1415926536E00
    C_SPEED_ANGSTROMS = 2.9979e18 # velocity of light in Ansgtroms/s
    C_SPEED_PI= 2.0 * np.pi * C_SPEED_ANGSTROMS

    PR0_exp= 0.0898 * (DENS**(1./6.))/np.sqrt(TEMP) #=(r0/debye)
    F00_exp=1.25E-9 * (DENS**(2./3.)) # normal field value in ues

    if PR0_exp > 1.:
        raise Exception('The plasma is too strongly correlated\ni.e. r0/debye=0.1\nthe line cannot be computed presently')

    # each input file contains tables for a given density
    # but 10 temperatures.
    # the parameter pr0 = r0/debye= 0.0898 Ne^{1/6} T^{-1/2},
    # (Ne in cm^-3, T in K) . For pr0 >1, we do not give a tabulation
    # yet now

    #id_max = 1

    ext = '.dat'

    for id in np.arange(id_max + 1):

        f = open(os.path.join(dir, 'profil' + str(id + 1) + ext))
        g = open(os.path.join(dir, 'index' + str(id + 1) + ext))

        # input tables Lyman/Balmer./Paschen, and indexes for reading inside
        # the tabulations which are given in dalfa units
        # dalfa=dlambda(Angstroms)/F00(ues)
        # corespond to delta alpha >0 (ie delta omega=domega <0)
        # Note that I(-domega) =I(domega
        # F00 = 1.25 10^{-9} Ne^{2/3}  (cm-3,ues) = normal field strength
        # fainu = wings factor in alfa units
        # ie (Idalfa)_wings= fainu/(dalfa**2.5)
        # (but not too far in the wings)

        # Regular expression for parsing the header
        st = '-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?'

        # Regular expression for parsing the tabulated raw_data
        st2 = '\s-?[0-9]+[0-9]*.?[0-9]*E-?\+?[0-9]+\s'

        tmp = f.readline()
        tmp2 = re.findall(st, tmp)
        N = np.int(tmp2[0])  # lower principal quantum number
        NP = np.int(tmp2[1])  # upper principal quantum number
        olam0 = np.float(tmp2[2])  # centre wavelength (A)

        tmp = f.readline()
        tmp2 = re.findall(st2,tmp)
        dense[id] = np.float(tmp2[0])
        #print 'dense[id]',dense[id]

        tmp = f.readline()
        tmp2 = re.findall(st2,tmp)
        f00[id] = np.float(tmp2[0])   # = 1.25 10^{-9} Ne^{2/3}

        tmp = f.readline()
        tmp2 = re.findall(st2,tmp)
        fainu = np.float(tmp2[0])

        ON = N
        ONP = NP
        DNU = 1.0 / ON **2 - 1.0/ ONP ** 2
        ambda = 911.7633455 * (ON * ONP) ** 2/ ((ONP - ON) * (ONP + ON))

        omega = C_SPEED_PI / ambda

        g.readline()
        g.readline()
        g.readline()
        g.readline()
        g.readline()

        if id < 17:
            #print 'id > 17 not triggered'
            f.readline()
            tmp = f.readline()
            tempe[0] = fort_conv(tmp[15:30])
            tempe[1] = fort_conv(tmp[30:55])
            tempe[2] = fort_conv(tmp[55:80])

            tmp = f.readline()
            pr0[id,0] = fort_conv(tmp[15:30])
            pr0[id,1] = fort_conv(tmp[30:55])
            pr0[id,2] = fort_conv(tmp[55:80])

            tmp = f.readline()
            dl12s[0] = fort_conv(tmp[15:30])
            dl12s[1] = fort_conv(tmp[30:55])
            dl12s[2] = fort_conv(tmp[55:80])

            tmp = f.readline()
            dl12[0] = fort_conv(tmp[15:30])
            dl12[1] = fort_conv(tmp[30:55])
            dl12[2] = fort_conv(tmp[55:80])
            f.readline()

            tmp = g.readline()
            ifm0 = np.int(tmp[9:12])
            ifm1 = np.int(tmp[13:16])

            itot=ifm1-ifm0+1
            #print 'ifm0,ifm1,itot',ifm0,ifm1,itot

            for i in np.arange(itot):
                tmp = f.readline()
                din[0,i] = fort_conv(tmp[0:11])
                sprof[0,i] = fort_conv(tmp[12:22])
                sprofs[0,i] = fort_conv(tmp[24:35])
                sprof[1,i] = fort_conv(tmp[37:47])
                sprofs[1,i] = fort_conv(tmp[49:59])
                sprof[2,i] = fort_conv(tmp[62:72])
                sprofs[2,i] = fort_conv(tmp[74:84])
                #print 'i,sprof 1, 2:',i,sprof[0,i],sprof[1,i]

            for jj in np.arange(max_d):
                din[1,jj] = din[0,jj]
                din[2,jj] = din[0,jj]

            # jtot(id,it)= number of wave lengths for the couple (T,Ne)
            for j in np.arange(3):
                tmp = False
                for i in np.arange(itot):
                    jtot[id,j] = i+1
                    if sprof[j,i] == 0.0:
                        jtot[id,j] = i#-1
                        tmp = True
                        break

        # Reference location 1717 on line 302
        f.readline()
        tmp = f.readline()
        tempe[3] = fort_conv(tmp[15:30])
        tempe[4] = fort_conv(tmp[30:55])
        tempe[5] = fort_conv(tmp[55:80])

        tmp = f.readline()
        pr0[id,3] = fort_conv(tmp[15:30])
        pr0[id,4] = fort_conv(tmp[30:55])
        pr0[id,5] = fort_conv(tmp[55:80])

        tmp = f.readline()
        dl12s[3] = fort_conv(tmp[15:30])
        dl12s[4] = fort_conv(tmp[30:55])
        dl12s[5] = fort_conv(tmp[55:80])

        tmp = f.readline()
        dl12[3] = fort_conv(tmp[15:30])
        dl12[4] = fort_conv(tmp[30:55])
        dl12[5] = fort_conv(tmp[55:80])
        f.readline()

        tmp = g.readline()
        ifm0 = np.int(tmp[9:12])
        ifm1 = np.int(tmp[13:16])

        itot=ifm1-ifm0+1
        #print 'ifm0,ifm1,itot',ifm0,ifm1,itot

        for i in np.arange(itot):
            tmp = f.readline()
            din[3,i] = fort_conv(tmp[0:11])
            sprof[3,i] = fort_conv(tmp[12:22])
            sprofs[3,i] = fort_conv(tmp[24:35])
            sprof[4,i] = fort_conv(tmp[37:47])
            sprofs[4,i] = fort_conv(tmp[49:59])
            sprof[5,i] = fort_conv(tmp[62:72])
            sprofs[5,i] = fort_conv(tmp[74:84])
            #print 'i,sprof 4, 5:',i,sprof[3,i],sprof[4,i]

        for jj in np.arange(max_d):
            din[4,jj] = din[3,jj]
            din[5,jj] = din[3,jj]

        for j in [3,4,5]:
            tmp = False
            for i in np.arange(itot):
                jtot[id,j] = i+1
                if sprof[j,i] == 0.0:
                    jtot[id,j] = i  #-1
                    tmp = True
                    break

        # Line 333
        f.readline()
        tmp = f.readline()
        tempe[6] = fort_conv(tmp[15:30])
        tempe[7] = fort_conv(tmp[30:55])
        tempe[8] = fort_conv(tmp[55:80])

        tmp = f.readline()
        pr0[id,6] = fort_conv(tmp[15:30])
        pr0[id,7] = fort_conv(tmp[30:55])
        pr0[id,8] = fort_conv(tmp[55:80])

        tmp = f.readline()
        dl12s[6] = fort_conv(tmp[15:30])
        dl12s[7] = fort_conv(tmp[30:55])
        dl12s[8] = fort_conv(tmp[55:80])

        tmp = f.readline()
        dl12[6] = fort_conv(tmp[15:30])
        dl12[7] = fort_conv(tmp[30:55])
        dl12[8] = fort_conv(tmp[55:80])
        f.readline()

        tmp = g.readline()
        ifm0 = np.int(tmp[9:12])
        ifm1 = np.int(tmp[13:16])

        itot=ifm1-ifm0+1
        #print 'ifm0,ifm1,itot',ifm0,ifm1,itot

        for i in np.arange(itot):
            tmp = f.readline()
            din[6,i] = fort_conv(tmp[0:11])
            sprof[6,i] = fort_conv(tmp[12:22])
            sprofs[6,i] = fort_conv(tmp[24:35])
            sprof[7,i] = fort_conv(tmp[37:47])
            sprofs[7,i] = fort_conv(tmp[49:59])
            sprof[8,i] = fort_conv(tmp[62:72])
            sprofs[8,i] = fort_conv(tmp[74:84])
            #print 'i,sprof 7, 8:',i,sprof[6,i],sprof[7,i]

        for jj in np.arange(max_d):
            din[7,jj] = din[6,jj]
            din[8,jj] = din[6,jj]

        for j in [6,7,8]:
            tmp = False
            for i in np.arange(itot):
                jtot[id,j] = i+1
                #print 'j,i,jtot(id,j)',j,i,jtot[id,j]
                if sprof[j,i] == 0.0:
                    jtot[id,j] = i#-1
                    tmp = True
                    break

        f.readline()
        tmp = f.readline()
        tempe[9] = fort_conv(tmp[15:30])

        tmp = f.readline()
        pr0[id,9] = fort_conv(tmp[15:30])

        tmp = f.readline()
        dl12s[9] = fort_conv(tmp[15:30])

        tmp = f.readline()
        dl12[9] = fort_conv(tmp[15:30])
        tmp = f.readline()

        tmp = g.readline()
        ifm0 = np.int(tmp[9:12])
        ifm1 = np.int(tmp[13:16])

        itot=ifm1-ifm0+1
        #print 'ifm0,ifm1,itot',ifm0,ifm1,itot

        for i in np.arange(itot):
            tmp = f.readline()
            din[9,i] = fort_conv(tmp[0:11])
            sprof[9,i] = fort_conv(tmp[12:22])
            sprofs[9,i] = fort_conv(tmp[24:35])
            #print 'i,sprof 10:',i,sprof[9,i]

        for j in [9]:
            tmp = False
            for i in np.arange(itot):
                jtot[id,j] = i+1
                if sprof[j,i] == 0.0:
                    jtot[id,j] = i#-1
                    tmp = True
                    break

        # Overwrite the temperatures with hard-coded values for some reason (line 388)
        tempe[0] =   2500.
        tempe[1] =   5000.
        tempe[2] =  10000.
        tempe[3] =  19950.
        tempe[4] =  39810.
        tempe[5] =  79430.
        tempe[6] = 158500.
        tempe[7] = 316200.
        tempe[8] = 631000.
        tempe[9] =1259600.

        for it in np.arange(10):
            if pr0[id,it] == 0.0: # non tabulated case
                pr0[id,it] = 0.0898*(dense[id]**(1./6.))/np.sqrt(tempe[it])

        PR00=0.0898*(dense[id]**(1./6.))/np.sqrt(TEMP)

        # conversion from alfa=dlambda(Angstroms)/F00 units
        # to normalized domega units (domega/F00) in rd/(s,ues),
        # normalized intensities are in units of (ues*s)/rd

        #jtot[id,:] = jtot[id,:]+1  # Added to agree with Fortran version

        for j in np.arange(10):
            for i in np.arange(int(jtot[id,j])):
                #print 'i,j,jtot[id,j]',i,j,jtot[id,j]
                # detunings for 1-np transition (alfa, omega, lambda units)
                dalfa=din[j,i]
                dlambda=f00[id]*dalfa   # angstroms
                domega= -C_SPEED_PI * dlambda / ((ambda + dlambda) * ambda) #rd/s
                dom[id,j,i]=domega/f00[id] # (rd/(s*ues)
                dom[id,j,i]=np.abs(dom[id,j,i])
                otrans= -C_SPEED_PI / (ambda * ambda)

                o1lines[id,j,i]=sprofs[j,i]
                o1line[id,j,i] = sprof[j,i]
                d1om[id,j,i]   =din[j,i]

                olines[id,j,i]=sprofs[j,i]/np.abs(otrans)
                oline[id,j,i] =sprof[j,i]/np.abs(otrans)

                #print '%.4e %.4e %.4e %.4e %.4e' %(tempe[j],dense[id],dom[id,j,i],olines[id,j,i],oline[id,j,i])

        #  asymptotic wing factor in normalized omega units
        # in the line wings : I(domega)=fainom/(domega**2.5)
        # In these units (rd/(s*ues)), the asymptotic constant
        # fainom is independnt of Ne and T
        # in the wings, one has I=fainom/(dom**2.5)
        fainom=fainu * (np.abs(otrans)) ** 1.5
        # in true angular frequency units one has
        # in the wings, one has I=fainom_exp/(domega**2.5)
        fainom_exp=fainom*(F00_exp**1.5)
        # in frequency units (1/s) nu=omega/(2*pi)
        # in the wings, one has I=fainum_exp/(dnu**2.5)
        fainum_exp=fainom_exp/( (np.pi * 2.)**1.5)

        #print 'fainum_exp',fainum_exp

        #print '-----'

        f.close()
        g.close()


    if write_netcdf:
        prefix = 'n_'+str(NP) + '_' + str(N) + '_'

        ncf.createDimension(prefix + 'id_maxi', id_maxi)

        wr_strs = ['N', 'NP', 'id_maxi', 'max_d', 'id_max',
                   'olam0', 'fainom', 'fainu', 'dense', 'tempe',
                   'f00', 'dl12', 'dl12s', 'pr0', 'jtot',
                   'dom', 'd1om', 'o1lines', 'o1line']

        wr_vars = [N, NP, id_maxi, max_d, id_max,
                   olam0, fainom, fainu, dense, tempe,
                   f00, dl12, dl12s, pr0, jtot,
                   dom, d1om, o1lines, o1line]

        wr_types = ['i', 'i', 'i', 'i', 'i',
                    'd', 'd', 'd', 'd', 'd',
                    'd', 'd', 'd', 'd', 'd',
                    'd', 'd', 'd', 'd']

        wr_dims = [('num',), ('num',), ('num',), ('num',), ('num',), ('num',), ('num',), ('num',), (prefix + 'id_maxi',), ('ten',), (prefix + 'id_maxi',), ('ten',), ('ten',), (prefix + 'id_maxi', 'ten',), (prefix + 'id_maxi', 'ten',), (prefix + 'id_maxi', 'ten', 'max_d',), (prefix + 'id_maxi', 'ten', 'max_d',), (prefix + 'id_maxi', 'ten', 'max_d',), (prefix + 'id_maxi', 'ten', 'max_d',)]

        for wr_str, wr_var, wr_type, wr_dim in zip(wr_strs, wr_vars, wr_types, wr_dims):
            variable = ncf.createVariable(prefix + wr_str, wr_type, wr_dim)
            variable[:] = wr_var

        # N_write = ncf.createVariable(prefix + 'N', 'i', ('num',))
        # N_write[:] = N
        # NP_write = ncf.createVariable(prefix + 'NP', 'd', ('num',))
        # NP_write[:] = NP
        # id_maxi_write = ncf.createVariable(prefix + 'id_maxi', 'd', ('num',))
        # id_maxi_write[:] = 30
        # max_d_write = ncf.createVariable(prefix+'max_d', 'd', ('num',))
        # max_d_write[:] = 60
        # id_max_write = ncf.createVariable(prefix+'id_max', 'd', ('num',))
        # id_max_write[:] = id_max
        # olam0_write = ncf.createVariable(prefix+'olam0', 'd', ('num',))
        # olam0_write[:] = olam0
        # fainom_write = ncf.createVariable(prefix+'fainom', 'd', ('num',))
        # fainom_write[:] = fainom
        # fainu_write = ncf.createVariable(prefix+'fainu', 'd', ('num',))
        # fainu_write[:] = fainu
        #
        # dense_write = ncf.createVariable(prefix+'dense','d',(prefix+'id_maxi',))
        # dense_write[:] = dense
        #
        # tempe_write = ncf.createVariable(prefix+'tempe','d',('ten',))
        # tempe_write[:] = tempe
        #
        # f00_write = ncf.createVariable(prefix+'f00','d',(prefix + 'id_maxi',))
        # f00_write[:] = f00
        #
        # dl12_write = ncf.createVariable(prefix+'dl12','d',('ten',))
        # dl12_write[:] = dl12
        #
        # dl12s_write = ncf.createVariable(prefix+'dl12s','d',('ten',))
        # dl12s_write[:] = dl12s
        #
        # pr0_write = ncf.createVariable(prefix+'pr0','d',(prefix+'id_maxi','ten',))
        # pr0_write[:] = pr0
        #
        # jtot_write = ncf.createVariable(prefix+'jtot','d',(prefix+'id_maxi','ten',))
        # jtot_write[:] = jtot
        #
        # dom_write = ncf.createVariable(prefix+'dom','d',(prefix+'id_maxi','ten','max_d',))
        # dom_write[:] = dom
        #
        # d1om_write = ncf.createVariable(prefix+'d1om','d',(prefix+'id_maxi','ten','max_d',))
        # d1om_write[:] = d1om
        #
        # o1lines_write = ncf.createVariable(prefix+'o1lines','d',(prefix+'id_maxi','ten','max_d',))
        # o1lines_write[:] = o1lines
        #
        # o1line_write = ncf.createVariable(prefix+'o1line','d',(prefix+'id_maxi','ten','max_d',))
        # o1line_write[:] = o1line



if write_netcdf:
    ncf.close()

    a = 5