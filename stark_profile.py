# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 12:11:41 2018

@author: jrh
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Feb 03 17:18:06 2018

@author: jrh
"""

#from scipy.interpolate import interp1d
#from scipy.optimize import curve_fit

import numpy as np
from scipy.constants import c, e, k
import time
import pystark




def data_tables(n_upper,n_lower):

    return tempe, olam0, N, NP, id_max, fainom, dense, f00, \
    dl12, dl12s, fainu, pr0, jtot, dom, d1om, o1line, o1lines


def fort_conv(string):
    if bool(string.strip()):
        return np.float(string)
    else:
        return 0.0


def piksrt(N,inarr):
    # Sorting by straight insertion, from Nummerical Recipes

    arr = np.zeros(N)
    arr[:] = inarr[:]
    
    for j in np.arange(1,N):
        a = arr[j]
        i = j-1
        while (i>0 and arr[i]>a):
            arr[i+1] = arr[i]
            i = i-1
        arr[i+1] = a
    return arr


def FINTRP(x1,x2,x3,y1,y2,y3,x):

#    return np.interp(x,[x1,x2,x3],[y1,y2,y3])

    if x == x2:
        return y2

    a12=x1-x2
    a22=x1-x3
    v1=y1-y2
    v2=y1-y3

    if ((y1 < y2) and (y2 < y3)) or ((y1 > y2) and (y2 > y3)):
        deter=v1*a22-v2*a12
        if np.abs(deter) < 1.0E-40:
            return y1+(x-x1)*(y3-y1)/(x3-x1)
        a21=x1*y1
        a11=a21-x2*y2
        a21=a21-x3*y3
        c=(a22*a11-a12*a21)/deter
        a=(-v2*a11+v1*a21)/deter
        b=(y1-a)*(x1-c)
        return a+b/(x-c)
    else:
        x1c=x1*x1
        a11=x1c-x2*x2
        a21=x1c-x3*x3
        deter=a11*a22-a12*a21
        if np.abs(deter) < 1.0E-40:
            raise Exception('FINTRP error: incorrect inputs')
        a=(a22*v1-a12*v2)/deter
        b=(-a21*v1+a11*v2)/deter
        c=y1-a*x1c-b*x1
        return (a*x+b)*x+c


def Voigt(x, amp, alpha, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    
    from scipy.special import wofz
    
    sigma = alpha / np.sqrt(2 * np.log(2))
    
    ans = np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)

    return amp*ans/np.max(ans)


def stark_profile(n_upper, n_lower, temperature, density):

    print('pystark')
    
    start = time.time()
    
    id_maxi = 30  # Maximum number of densities
    max_d = 60  # Maximum number of detunings
    
    TEMP = temperature * e / k  # temperature in K
    DENS = density * 1.0e-6  # electronic density in cm-3
    
    # ------------------------------------------------------------------------
    # Data tables
    # ------------------------------------------------------------------------
    if n_lower == 1:
        # Lyman series
        if n_upper == 2:
            from pystark.stark_2_1_data import stark_2_1_data as sdata
        if n_upper == 3:
            from pystark.stark_3_1_data import stark_3_1_data as sdata
        if n_upper == 4:
            from pystark.stark_4_1_data import stark_4_1_data as sdata
        if n_upper == 5:
            from pystark.stark_5_1_data import stark_5_1_data as sdata
        if n_upper == 6:
            from pystark.stark_6_1_data import stark_6_1_data as sdata
        if n_upper == 7:
            from pystark.stark_7_1_data import stark_7_1_data as sdata
        if n_upper == 8:
            from pystark.stark_8_1_data import stark_8_1_data as sdata
        if n_upper == 9:
            from pystark.stark_9_1_data import stark_9_1_data as sdata
        if n_upper == 10:
            from pystark.stark_10_1_data import stark_10_1_data as sdata
        if n_upper == 11:
            from pystark.stark_11_1_data import stark_11_1_data as sdata
        if n_upper == 12:
            from pystark.stark_12_1_data import stark_12_1_data as sdata
        if n_upper == 13:
            from pystark.stark_13_1_data import stark_13_1_data as sdata
        if n_upper == 14:
            from pystark.stark_14_1_data import stark_14_1_data as sdata
        if n_upper == 15:
            from pystark.stark_15_1_data import stark_15_1_data as sdata
        if n_upper == 16:
            from pystark.stark_16_1_data import stark_16_1_data as sdata
        if n_upper == 17:
            from pystark.stark_17_1_data import stark_17_1_data as sdata
        if n_upper == 18:
            from pystark.stark_18_1_data import stark_18_1_data as sdata
        if n_upper == 19:
            from pystark.stark_19_1_data import stark_19_1_data as sdata
        if n_upper == 20:
            from pystark.stark_20_1_data import stark_20_1_data as sdata
        if n_upper == 21:
            from pystark.stark_21_1_data import stark_21_1_data as sdata
        if n_upper == 22:
            from pystark.stark_22_1_data import stark_22_1_data as sdata
        if n_upper == 23:
            from pystark.stark_23_1_data import stark_23_1_data as sdata
        if n_upper == 24:
            from pystark.stark_24_1_data import stark_24_1_data as sdata
        if n_upper == 25:
            from pystark.stark_25_1_data import stark_25_1_data as sdata
        if n_upper == 26:
            from pystark.stark_26_1_data import stark_26_1_data as sdata
        if n_upper == 27:
            from pystark.stark_27_1_data import stark_27_1_data as sdata
        if n_upper == 28:
            from pystark.stark_28_1_data import stark_28_1_data as sdata
        if n_upper == 29:
            from pystark.stark_29_1_data import stark_29_1_data as sdata
        if n_upper == 30:
            from pystark.stark_30_1_data import stark_30_1_data as sdata
    if n_lower == 2:
        # Balmer series
        if n_upper == 3:    
            from pystark.stark_3_2_data import stark_3_2_data as sdata
        if n_upper == 4:
            from pystark.stark_4_2_data import stark_4_2_data as sdata
        if n_upper == 5:
            from pystark.stark_5_2_data import stark_5_2_data as sdata
        if n_upper == 6:
            from pystark.stark_6_2_data import stark_6_2_data as sdata
        if n_upper == 7:
            from pystark.stark_7_2_data import stark_7_2_data as sdata
        if n_upper == 8:
            from pystark.stark_8_2_data import stark_8_2_data as sdata
        if n_upper == 9:
            from pystark.stark_9_2_data import stark_9_2_data as sdata
        if n_upper == 10:
            from pystark.stark_10_2_data import stark_10_2_data as sdata
        if n_upper == 11:
            from pystark.stark_11_2_data import stark_11_2_data as sdata
        if n_upper == 12:
            from pystark.stark_12_2_data import stark_12_2_data as sdata
        if n_upper == 13:
            from pystark.stark_13_2_data import stark_13_2_data as sdata
        if n_upper == 14:
            from pystark.stark_14_2_data import stark_14_2_data as sdata
        if n_upper == 15:
            from pystark.stark_15_2_data import stark_15_2_data as sdata
        if n_upper == 16:
            from pystark.stark_16_2_data import stark_16_2_data as sdata
        if n_upper == 17:
            from pystark.stark_17_2_data import stark_17_2_data as sdata
        if n_upper == 18:
            from pystark.stark_18_2_data import stark_18_2_data as sdata
        if n_upper == 19:
            from pystark.stark_19_2_data import stark_19_2_data as sdata
        if n_upper == 20:
            from pystark.stark_20_2_data import stark_20_2_data as sdata
        if n_upper == 21:
            from pystark.stark_21_2_data import stark_21_2_data as sdata
        if n_upper == 22:
            from pystark.stark_22_2_data import stark_22_2_data as sdata
        if n_upper == 23:
            from pystark.stark_23_2_data import stark_23_2_data as sdata
        if n_upper == 24:
            from pystark.stark_24_2_data import stark_24_2_data as sdata
        if n_upper == 25:
            from pystark.stark_25_2_data import stark_25_2_data as sdata
        if n_upper == 26:
            from pystark.stark_26_2_data import stark_26_2_data as sdata
        if n_upper == 27:
            from pystark.stark_27_2_data import stark_27_2_data as sdata
        if n_upper == 28:
            from pystark.stark_28_2_data import stark_28_2_data as sdata
        if n_upper == 29:
            from pystark.stark_29_2_data import stark_29_2_data as sdata
        if n_upper == 30:
            from pystark.stark_30_2_data import stark_30_2_data as sdata
    if n_lower == 3:
        # Paschen series
        if n_upper == 4:
            from pystark.stark_4_3_data import stark_4_3_data as sdata
        if n_upper == 5:
            from pystark.stark_5_3_data import stark_5_3_data as sdata
        if n_upper == 6:
            from pystark.stark_6_3_data import stark_6_3_data as sdata
        if n_upper == 7:
            from pystark.stark_7_3_data import stark_7_3_data as sdata
        if n_upper == 8:
            from pystark.stark_8_3_data import stark_8_3_data as sdata
        if n_upper == 9:
            from pystark.stark_9_3_data import stark_9_3_data as sdata
        if n_upper == 10:
            from pystark.stark_10_3_data import stark_10_3_data as sdata
        if n_upper == 11:
            from pystark.stark_11_3_data import stark_11_3_data as sdata
        if n_upper == 12:
            from pystark.stark_12_3_data import stark_12_3_data as sdata
        if n_upper == 13:
            from pystark.stark_13_3_data import stark_13_3_data as sdata
        if n_upper == 14:
            from pystark.stark_14_3_data import stark_14_3_data as sdata
        if n_upper == 15:
            from pystark.stark_15_3_data import stark_15_3_data as sdata
        if n_upper == 16:
            from pystark.stark_16_3_data import stark_16_3_data as sdata
        if n_upper == 17:
            from pystark.stark_17_3_data import stark_17_3_data as sdata
        if n_upper == 18:
            from pystark.stark_18_3_data import stark_18_3_data as sdata
        if n_upper == 19:
            from pystark.stark_19_3_data import stark_19_3_data as sdata
        if n_upper == 20:
            from pystark.stark_20_3_data import stark_20_3_data as sdata
        if n_upper == 21:
            from pystark.stark_21_3_data import stark_21_3_data as sdata
        if n_upper == 22:
            from pystark.stark_22_3_data import stark_22_3_data as sdata
        if n_upper == 23:
            from pystark.stark_23_3_data import stark_23_3_data as sdata
        if n_upper == 24:
            from pystark.stark_24_3_data import stark_24_3_data as sdata
        if n_upper == 25:
            from pystark.stark_25_3_data import stark_25_3_data as sdata
        if n_upper == 26:
            from pystark.stark_26_3_data import stark_26_3_data as sdata
        if n_upper == 27:
            from pystark.stark_27_3_data import stark_27_3_data as sdata
        if n_upper == 28:
            from pystark.stark_28_3_data import stark_28_3_data as sdata
        if n_upper == 29:
            from pystark.stark_29_3_data import stark_29_3_data as sdata
        if n_upper == 30:
            from pystark.stark_30_3_data import stark_30_3_data as sdata

    
    tempe, olam0, N, NP, id_max, fainom, dense, f00, dl12, dl12s, fainu, pr0, \
    jtot, dom, d1om, o1line, o1lines = sdata()

    load_time = time.time() - start
    print('load_time:', load_time)

    # tempe     - electron temperature (K)
    # olam0     - line centre wavlength (A)
    # N         - lower principal quantum number
    # NP        - upper principal quantum number
    # id_max    -
    # fainom    -
    # dense     - electron density (cm ** -3)
    # f00       - normal Holtsmark field strength (30 kV / m)
    # dl12      -
    # dl12s     -
    # fainu     - Asymptotic value of iStark * (alpha ** 2.5)
    # pr0       - Ratio of the mean interelectronic distance to the electronic Debye length
    # jtot      -
    # dom       -
    # d1om      -
    # o1line    -
    # o1lines   -

    jtot = jtot.astype(np.int)
    
    #dom = np.zeros([id_maxi,10,max_d])
    oline = np.zeros([id_maxi, 10, max_d])
    olines = np.zeros([id_maxi, 10, max_d])
    
    #d1om = np.zeros([id_maxi,10,max_d])
    #o1line = np.zeros([id_maxi,10,max_d])
    #o1lines = np.zeros([id_maxi,10,max_d])
    
    domm = np.zeros(100000)
    dom0 = np.zeros(10000)
    tprof = np.zeros([id_maxi, 10, 10000])
    tprofs = np.zeros([id_maxi, 10, 10000])
    uprof = np.zeros([id_maxi, 10000])
    uprofs = np.zeros([id_maxi, 10000])

    cspeed = c * 1e10  # velocity of light in Ansgtroms/s
    cspeed_pi = 2.0 * np.pi * cspeed
    
    PR0_exp = 0.0898 * (DENS ** (1./6.))/ np.sqrt(TEMP) #=(r0/debye)
    F00_exp = 1.25E-9 * (DENS ** (2./3.))  # normal field value in ues
    
    ON=N
    ONP=NP
    #DNU=1.0/ON**2 - 1.0/ONP**2
    ambda=911.7633455*(ON*ONP)**2/((ONP-ON)*(ONP+ON))
    
    omega=cspeed_pi/ambda
    otrans=-cspeed_pi/(ambda*ambda)
    
    olines=o1lines/np.abs(otrans)
    oline=o1line/np.abs(otrans)
    
    if PR0_exp > 1.:
        raise Exception('The plasma is too strongly correlated\ni.e. r0/debye=0.1\nthe line cannot be computed presently')
    
    #fainom_exp=fainom*(F00_exp**1.5)
    #fainum_exp=fainom_exp/( (OPI*2.)**1.5)
    
    #========================
    # TABULATION Format CDS
    #   si on veut ecrire 
    #  n -np lambda0 kalpha Ne E0 T R0/Debye Dalpha iDoppler iStark
    
    #IN_cds= N+0.01
    #INP_cds = NP+0.01
    
    #***********************************************************
    # Don't edit the CDS format...
    #***********************************************************
    
    # Skipped the code in the IF statement starting at line 470, since it
    # isn't used, if (.FALSE.) ...
    
    #========================================================
    # change sligtly the value of the input density
    # DENS in order to remain , as far as possible, inside the tabulation
    
    if np.abs(DENS-dense[0])/DENS <= 1.0E-3:
        DENS = dense[0]*1.001
        
    for id in np.arange(1,id_max+1):
        if np.abs(DENS-dense[id])/DENS <= 1.0E-3:
            DENS = dense[id]*0.999
            
    #==============================================
    # define an unique detunings grid - domm -  for the tabulated
    # profiles ( various temperatures , densities)
    # calculate all the line shapes for this  common grid
    # units used at this points are Domega_new= Delta(omega)/F00
    #                                      in rd/(s-1 ues)
            
    inc = 0
    domm[inc] = 0.0
    # ---- Look to replace this loop
    for id in np.arange(id_max+1):
        for j in np.arange(10):
            #print 'jtot[id,j]',jtot[id,j]
            for i in np.arange(1, jtot[id, j]):
                inc=inc+1
                dom0[inc] = dom[id,j,i]
    
    inc = np.count_nonzero(dom)
    npik = inc + 1
    #nut=10000
    
    # Calling numpy sort instead of piksrt
    tmp = np.sort(dom0[0:npik])
    dom0[0:npik] = tmp[0:npik]
    # dom0 seems to agree with the FORTRAN version
    
    inc = 0
    domm[0] = 0.0
    #print 'npik',npik
    # ---- Look to replace this loop
    for i in np.arange(1,npik):
        #print 'i',i
        dif = (dom0[i]-dom0[i-1])
        if dif <= 1.0E-6:
            continue
        if dif/np.abs(dom0[i]) <= 0.1:
            continue
        inc = inc+1
        domm[inc] = dom0[i]
        
    jdom = inc + 1  # One line after marker 35
    
    for id in np.arange(id_max):
        for j in np.arange(10):
            if pr0[id,j] > 1.0:
                continue
            
            tprof[id,j,0] = oline[id,j,0]
            tprofs[id,j,0] = olines[id,j,0]
            
            #print 'id,j',id,j
            #print 'tprof,tprofs:',tprof[id,j,0],tprofs[id,j,0]
    
            if jtot[id,j] == 0:
                continue
            
            for i in np.arange(1,jdom+1):
                skip1 = False
                skip2 = False
                #print 'i',i
                domeg=domm[i]
                ij_max=jtot[id,j]
                #print 'domeg,ij_max',domeg,ij_max
                for ij in np.arange(1,ij_max-1):
                    #print 'ij',ij
                    test=(domeg-dom[id,j,ij])*(domeg-dom[id,j,ij-1])
                    #print 'test1:',test
                    if test <= 0.0:
                        #print 'triggered test1'
                        x1=dom[id,j,ij-1]
                        x2=dom[id,j,ij]
                        x3=dom[id,j,ij+1]
                        y1=oline[id,j,ij-1]
                        y2=oline[id,j,ij]
                        y3=oline[id,j,ij+1]
                        #print 'x1,x2,x3',x1,x2,x3
                        #print 'y1,y2,y3',y1,y2,y3
                        tprof[id,j,i] = FINTRP(x1,x2,x3,y1,y2,y3,domeg)
                        y1=olines[id,j,ij-1]
                        y2=olines[id,j,ij]
                        y3=olines[id,j,ij+1]
                        tprofs[id,j,i] = FINTRP(x1,x2,x3,y1,y2,y3,domeg)
                        #print 'tprof[id,j,i]',tprof[id,j,i]
                        #print 'tprofs[id,j,i]',tprofs[id,j,i]
                        skip1 = True
                        skip2 = True
                        break
                    
                if skip1 is False:
                    test=(domeg-dom[id,j,ij_max-2])*(domeg-dom[id,j,ij_max-1])
                    #print 'test2:',test
                    #print 'domeg',domeg
                    #print 'dom[id,j,ij_max-1]',dom[id,j,ij_max-2]
                    #print 'dom[id,j,ij_max]',dom[id,j,ij_max-1]
                    if test <= 0.0:
                        #print 'triggered test2'
                        x1=dom[id,j,ij_max-3]
                        x2=dom[id,j,ij_max-2]
                        x3=dom[id,j,ij_max-1]
                        y1=oline[id,j,ij_max-3]
                        y2=oline[id,j,ij_max-2]
                        y3=oline[id,j,ij_max-1]
                        tprof[id,j,i]= FINTRP(x1,x2,x3,y1,y2,y3,domeg)
                        y1=olines[id,j,ij_max-3]
                        y2=olines[id,j,ij_max-2]
                        y3=olines[id,j,ij_max-1]
                        tprofs[id,j,i]= FINTRP(x1,x2,x3,y1,y2,y3,domeg)
                        skip2 = True
                        #print 'x1,x2,x3',x1,x2,x3
                        #print 'y1,y2,y3',y1,y2,y3
                        #print 'tprof[id,j,i]',tprof[id,j,i]
                        #print 'tprofs[id,j,i]',tprofs[id,j,i]
                        continue
                    
                if skip2 is False:
                    if domeg > dom[id,j,ij_max]:
                        #print 'triggered test3'
                        tprof[id,j,i] =fainom/(domeg**2.5)
                        tprofs[id,j,i]=tprof[id,j,i]
                        continue
    
    # We can skip writing the intermediate file
    
    if DENS >= 2.0*dense[id_max]:
        raise Exception('Your input density is higher than the largest tabulated value %f' %dense[id_max])
    
    if DENS <= dense[0]:
        raise Exception('Your input density is smaller than the smallest tabulated value %f' %dense[0])
        
    if TEMP >= tempe[9]:
        raise Exception('Your input temperature is higher than the largest tabulated value %f' %tempe[9])
        
    if TEMP <= tempe[0]:
        raise Exception('Your input temperature is lower than the smallest tabulated value %f' %tempe[0])
        
    
    for id in np.arange(id_max):
        otest_dens = (DENS-dense[id])*(DENS-dense[id+1])
        if otest_dens <= 0.0:
            dense1 = dense[id]
            dense2 = dense[id+1]
            id1 =id
            id2 = id+1
            break
    
    if DENS >= dense[id_max]:
        dense1 = dense[id_max-1]
        dense2 = dense[id_max]
        id1 = id_max-1
        id2 = id_max
        
    for it in np.arange(10):
        otest = (TEMP-tempe[it])*(TEMP-tempe[it+1])
        if otest <= 0.0:
            it1 = it
            it2 = it+1
            #pr01 = pr0[id2,it1] # max value of pr0 for T1,T2,dense1,dense2
            tempe1 = tempe[it]
            tempe2 = tempe[it+1]
            break
        
    #print 'id1,id2',id1,id2
    #print 'jdom:',jdom
    
    # interpolation in temperature          
    for id in np.arange(id1,id2+1):
        for i in np.arange(jdom):
            uprof[id,i] = tprof[id,it1,i] + (TEMP-tempe1) * (tprof[id,it2,i]-tprof[id,it1,i]) / (tempe2-tempe1)
            uprofs[id,i] = tprofs[id,it1,i] + (TEMP-tempe1) * (tprofs[id,it2,i]-tprofs[id,it1,i]) / (tempe2-tempe1)
            #print 'id,i:',id,i
            #print 'tprof[id,it,i]:',tprof[id,it,i]
            #print 'uprof[id,i]:',uprof[id,i]
            #print 'uprofs[id,i]:',uprofs[id,i]
            
#    print ' Interpolated Profiles I(delta_nu) in s'
#    print ' for detunings delta_nu in s-1' 
#    print '                 I(delta_nu)=I(-delta_nu)'
#    print ' first column is the corresponding detuning in Angstroms'
#    print ' second column is the detuning delta_nu'
#    print ' column 3: Doppler convolved profile in s'
#    print '                 I(delta_nu)=I(-delta_nu)'
#    print ' column 4:  pure Stark profile'
#    print ''
#    print '  Wings factor Knu in s**1.5 such that'
#    print '    I(delta_nu)=Knu/(delta_nu^2.5)'
#    print '      Knu=%g' %fainum_exp
    
    #print jdom
    
    delta_lambda = np.zeros(jdom)
    delta_nu = np.zeros(jdom)
    wprof_nu = np.zeros(jdom)
    wprofs_nu = np.zeros(jdom)
    
    for i in np.arange(jdom):
        wprof=uprof[id1,i]+(DENS-dense1)*(uprof[id2,i]-uprof[id1,i])/(dense2-dense1)
        wprofs=uprofs[id1,i]+(DENS-dense1)*(uprofs[id2,i]-uprofs[id1,i])/(dense2-dense1)
        delta_omega=domm[i]*F00_exp
        delta_nu[i]=delta_omega/(2 * np.pi)
        delta_lambda[i]= ambda*delta_omega/(omega + delta_omega)
        wprof_nu[i]=(wprof/F00_exp)*(2. * np.pi)
        wprofs_nu[i]=(wprofs/F00_exp)*(2. * np.pi)
#        print '%e %e %e %e' %(delta_lambda[i],delta_nu[i],wprof_nu[i],wprofs_nu[i])
        
    delta_lambda2 = np.concatenate((-delta_lambda[::-1],delta_lambda)) + olam0
#    delta_nu2 = np.concatenate((-delta_nu[::-1],delta_nu))
    wprof_nu2 = np.concatenate((wprof_nu[::-1],wprof_nu))
    wprofs_nu2 = np.concatenate((wprofs_nu[::-1],wprofs_nu))
    
    end = time.time()
    
    print('Time elapsed (s): ',end-start)
    
    return delta_lambda2, wprof_nu2, wprofs_nu2


if __name__ == '__main__':        
    
    x, y, ys = stark_profile(2, 1, 8.621738, 1.0E20)
    
    olam0 = np.mean(x)
    
    # Read output from the fortran code
    #f77res = np.loadtxt('C:/Temp/98A/ba/ba08/fort.74',skiprows=5)
    # f77res = np.loadtxt('C:/Temp/98A/pa/pa04/fort.74',skiprows=5)
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(x, y, '.-b')
    # plt.plot(olam0+f77res[:,0],f77res[:,2],'x-r')
    # plt.legend(['PyStark','Fortran'])
    
    #fit = curve_fit(Voigt,delta_lambda2,wprof_nu2)

    #plt.plot(delta_lambda2,Voigt(delta_lambda2,2.52E-11,0.235,0.03),'g')
    #xdat = np.linspace(-2.0,2.0,150)
    #plt.plot(xdat,Voigt(xdat,fit[0][0],fit[0][1],fit[0][2]),'g')
    plt.xlim(olam0 + np.array([-30.0, 30.0]))
    
    #diff = 100.0*np.abs(y-f77res[:,2])/f77res[:,2]
    #plt.figure()
    #plt.plot(f77res[:,1],diff)
    #plt.xlim([0.0,1.0])
    plt.xlabel('wavelength (A)')
    plt.show()