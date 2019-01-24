import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import wofz
from scipy.constants import physical_constants, atomic_mass, c, h, e, k, m_e
import scipy.signal
import pystark


class BalmerLineshape(object):
    def __init__(self, n_upper, dens, temp, bfield, viewangle=0., line_model='rosato', wl_axis=None,
                 wl_centre=None, npts=None, isotope='D', override_input_check=False):

        """ Hydrogen Balmer series spectral lineshape. Area-normalised to 1.
        
        Input / output is in wavelength space, while internal calculations are done in frequency space where possible /
        appropriate. Only the Balmer series is currently supported, however 'stehle' and 'stehle param' models do
        support some Paschen lines so it would be easy to incorporate those should you need to. Not optimised for
        speed yet.
        
        :param n_upper: upper principal quantum number of transition.
        :type n_upper: int

        :param dens: density [ m^-3 ]
        :type dens: float

        :param temp: temperature [ eV ]
        :type temp: float

        :param bfield: magnetic field strength [ T ]
        :type dens: float

        :param viewangle: [ rad ]
        :type viewangle: float

        :param line_model: can be:

        - 'voigt':

        - 'rosato':

        - 'stehle':

        - 'stehle_param':

        :type line_model: str

        :param wl_axis: [ m ]
        :type wl_axis: np.array

        :param wl_centre: [ m ]
        :type wl_centre: float

        :param npts: number of points in auto-generated wl_axis (if wl_axis not specified)
        :type npts: int

        :param isotope: can be 'H', 'D', 'T'
        :type isotope: str

        :param override_input_check: skip checking that the inputs lie within valid ranges - use with care
        :type override_input_check: bool
        """

        if npts is None:
            npts = 3001

        # if no wavelength centre supplied, retrieve NIST value
        if wl_centre is None:
            wl_centre = pystark.get_wl_centre(n_upper)

        # if no wavelength axis supplied, generate reasonable axis
        if wl_axis is None:
            wl_axis = pystark.get_wl_axis(n_upper, dens, temp, bfield, npts=npts, wl_centre=wl_centre)

        self.npts = npts
        self.line_model = line_model
        self.n_upper = n_upper
        self.n_lower = 2
        self.dens = dens
        self.temp = temp
        self.bfield = bfield
        self.viewangle = viewangle
        self.isotope = isotope
        self.mass = pystark.get_h_isotope_mass(isotope)
        self.wl_centre = wl_centre
        self.wl_axis = wl_axis

        if override_input_check is False:
            self.input_check()

        # frequency axis for internal use only
        self.freq_axis = pystark.get_freq_axis(n_upper, dens, temp, bfield, no_fwhm=30, npts=self.npts,
                                               wl_centre=wl_centre)

        self.freq_axis_conv = pystark.get_freq_axis_conv(self.freq_axis)
        self.freq_centre = c / wl_centre
        self.energy_centre = h * self.freq_centre  # [ eV ]

        # generate lineshape
        if self.line_model == 'rosato':
            ls_sz = self.make_rosato()

            # Calculate Doppler lineshape
            ls_d = pystark.doppler_lineshape(self.freq_axis_conv, self.freq_centre, self.temp, self.mass, x_units='Hz')

            # convolution in frequency space
            ls_szd = scipy.signal.fftconvolve(ls_sz, ls_d, 'same')  # [ / Hz ]
            ls_szd /= np.trapz(ls_szd, self.freq_axis)
        else:

            if self.line_model == 'voigt':
                ls_szd = self.make_voigt()

            elif self.line_model == 'stehle_param':
                ls_s = self.make_stehle_param()

                # Calculate Doppler lineshape
                ls_d = pystark.doppler_lineshape(self.freq_axis_conv, self.freq_centre, self.temp, self.mass, x_units='Hz')

                # convolution in frequency space
                ls_sd = scipy.signal.fftconvolve(ls_s, ls_d, 'same')  # [ / Hz ]
                ls_sd /= np.trapz(ls_sd, self.freq_axis)

                ls_szd = self.zeeman_split(self.freq_axis, self.freq_centre, ls_sd)

            elif self.line_model == 'stehle':
                ls_s = self.make_stehle()

                # Calculate Doppler lineshape
                ls_d = pystark.doppler_lineshape(self.freq_axis_conv, self.freq_centre, self.temp, self.mass, x_units='Hz')

                # convolution in frequency space
                ls_sd = scipy.signal.fftconvolve(ls_s, ls_d, 'same')  # [ / Hz ]
                ls_sd /= np.trapz(ls_sd, self.freq_axis)

                ls_szd = self.zeeman_split(self.freq_axis, self.freq_centre, ls_sd)

        # interpolate onto wavelength grid
        x_out, x_centre_out, self.ls_szd = pystark.convert_ls_units(self.freq_axis, self.freq_centre, mode='interp',
                                                                    x_out=self.wl_axis, ls=ls_szd)

        self.ls_szd = pystark.ls_norm(self.ls_szd, self.wl_axis)

    def input_check(self):
        """
        Check __init__ arguments are valid.

        """

        # TODO this all seems a bit convoluted...

        line_models = ['voigt', 'rosato', 'stehle', 'stehle_param',]
        n_upper_range = [None, (3, 7), (3, 30), (3, 9)]
        dens_range = [None, (1e19, 1e21), (1e16, 1e25), (0., 1e22)]
        temp_range = [None, (0.32, 32), (0.22, 110), (0., 1000)]
        bfield_range = [None, (0, 5), (0, 5), (0, 5)]

        param_ranges = list(zip(line_models, n_upper_range, dens_range, temp_range, bfield_range))
        columns = ['line model names', 'n upper range', 'dens range', 'temp range', 'bfield range']
        param_ranges = pd.DataFrame(data=param_ranges, columns=columns)

        assert sum(param_ranges['line model names'].isin([self.line_model]))

        n_upper_range = param_ranges['n upper range'][param_ranges['line model names'] == self.line_model].values[0]
        dens_range = param_ranges['dens range'][param_ranges['line model names'] == self.line_model].values[0]
        temp_range = param_ranges['temp range'][param_ranges['line model names'] == self.line_model].values[0]
        bfield_range = param_ranges['bfield range'][param_ranges['line model names'] == self.line_model].values[0]

        if n_upper_range is not None:
            assert (self.n_upper in range(n_upper_range[0], n_upper_range[1]))
        if dens_range is not None:
            assert (dens_range[0] <= self.dens <= dens_range[1])
        if temp_range is not None:
            assert (temp_range[0] <= self.temp <= temp_range[1])
        if bfield_range is not None:
            assert (bfield_range[0] <= self.bfield <= bfield_range[1])

    def make_rosato(self):
        """
        Stark-Zeeman lineshape interpolated using the Rosato et al. tabulated_data.
        tables.
        
        outputs on a uniform frequency grid.

        :return:

        """

        # use fwhm estimate to generate inputs to interpolation routine
        fwhm_estimate = pystark.estimate_fwhm(self.n_upper, self.dens, self.temp, self.bfield, isotope=self.isotope)
        no_fwhm = 20
        wmax = no_fwhm * fwhm_estimate * h / e  # [ eV ]

        # load interpolated Stark-Zeeman lineshape using python wrapper
        detunings_rosato, ls_sz = pystark.rosato_wrapper(self.n_upper, self.dens, self.temp, self.bfield,
                                                         self.viewangle, wmax, self.npts, display=False)

        freqs_rosato = e * detunings_rosato / h + self.freq_centre  # [ Hz ]
        ls_sz *= h / e  # [ / Hz ]

        # area normalise
        ls_sz /= np.trapz(ls_sz, freqs_rosato)

        # interpolate onto class frequency grid
        ls_sz = np.interp(self.freq_axis, freqs_rosato, ls_sz)

        return ls_sz

    def make_voigt(self):
        """
        Voigt profile is related to the real part of the Faddeeva function -- avoids numerical convolution of
        Lorentzian and Gaussian profiles.
        
        :return:

        """

        hwhm_doppler = pystark.doppler_fwhm(self.n_upper, self.temp, self.mass) / 2
        sigma_doppler = hwhm_doppler / np.sqrt(2 * np.log(2))

        hwhm_stark = pystark.stark_fwhm(self.n_upper, self.dens) / 2

        ls_sd = np.real(wofz(((self.freq_axis - self.freq_centre) + 1j * hwhm_stark) / sigma_doppler / np.sqrt(2))) / \
                sigma_doppler / np.sqrt(2 * np.pi)

        ls_szd = self.zeeman_split(self.freq_axis, self.freq_centre, ls_sd)

        return ls_szd

    def make_stehle_param(self):

        # Paramaterised MMM Stark profile coefficients from Bart's paper

        loman_ij_abc = {
            '32':  [0.7665, 0.064, 3.710e-18],
            '42':  [0.7803, 0.050, 8.425e-18],
            '52':  [0.6796, 0.030, 1.310e-15],
            '62':  [0.7149, 0.028, 3.954e-16],
            '72':  [0.7120, 0.029, 6.258e-16],
            '82':  [0.7159, 0.032, 7.378e-16],
            '92':  [0.7177, 0.033, 8.947e-16],
            '43':  [0.7449, 0.045, 1.330e-16],
            '53':  [0.7356, 0.044, 6.640e-16],
            '63':  [0.7118, 0.016, 2.481e-15],
            '73':  [0.7137, 0.029, 3.270e-15],
            '83':  [0.7133, 0.032, 4.343e-15],
            '93':  [0.7165, 0.033, 5.588e-15],
        }

        loman_abc = loman_ij_abc[str(self.n_upper) + str(self.n_lower)]
        a_ij, b_ij, c_ij = loman_abc

        delta_lambda_12ij = c_ij * (self.dens ** a_ij) / (self.temp ** b_ij)  # nm

        ls_s = 1 / (abs((self.wl_axis - self.wl_centre) * 1e9) ** (5. / 2.) +
                    (delta_lambda_12ij / 2) ** (5. / 2.))

        ls_s /= np.trapz(ls_s, self.wl_axis)

        x_out, x_centre_out, ls_s = pystark.convert_ls_units(self.wl_axis, self.wl_centre, mode='interp', x_out=self.freq_axis, ls=ls_s)

        return ls_s

    def make_stehle(self):
        """
        ensure given n_upper + n_lower fall within tabulated values

        :return:

        """

        temp_k = self.temp * e / k  # temperature in K
        dens_cm = self.dens * 1.e-6  # electronic density in cm-3
        prefix = 'n_' + str(self.n_upper) + '_' + str(self.n_lower) + '_'

        # extract raw tabulated tabulated_data
        tab_temp_k = np.array(pystark.nc.variables[prefix + 'tempe'].data)  # tabulated electron temperatures (K)
        olam0 = pystark.nc.variables[prefix + 'olam0'].data  # line centre wavelength (A)
        num_tab_dens = pystark.nc.variables[prefix + 'id_max'].data
        fainom = pystark.nc.variables[prefix + 'fainom'].data
        tab_dens_cm = np.array(pystark.nc.variables[prefix + 'dense'].data)  # tabulated electron densities  (cm ** -3)
        f00 = np.array(pystark.nc.variables[prefix + 'f00'].data)  # normal Holtsmark field strength (30 kV / m)
        dl12 = np.array(pystark.nc.variables[prefix + 'dl12'].data)
        dl12s = np.array(pystark.nc.variables[prefix + 'dl12s'].data)
        fainu = pystark.nc.variables[
            prefix + 'fainu'].data  # Asymptotic value of iStark * (alpha ** 2.5) ("wings factor in alfa units")
        pr0 = np.array(pystark.nc.variables[
                           prefix + 'pr0'].data)  # Ratio of the mean interelectronic distance to the electronic Debye length
        jtot = np.array(pystark.nc.variables[prefix + 'jtot'].data,
                        dtype=np.int)  # "number of wave lengths for the couple (T,Ne)"
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
            raise Exception(
                'Your input density is higher than the largest tabulated value %f' % tab_dens_cm[num_tab_dens])

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
        # rydberg_m = Rydberg / (1. + (electron_mass / physical_constants['deuteron mass'][0]))
        # wl_0_angst = 1e10 * (rydberg_m * (1 / n_lower ** 2 - 1 / n_upper ** 2)) ** -1

        wl_centre_angst = self.wl_centre * 1e10

        c_angst = c * 1e10  # velocity of light in Ansgtroms / s
        angular_freq_0 = 2 * np.pi * c_angst / wl_centre_angst  # rad / s

        otrans = -2 * np.pi * c_angst / wl_centre_angst ** 2

        olines = o1lines / np.abs(otrans)
        oline = o1line / np.abs(otrans)

        # Limit analysis_tools to uncorrelated plasmas.
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
                            tprof[id, j, i] = pystark.FINTRP(x1, x2, x3, y1, y2, y3, domeg)
                            y1 = olines[id, j, ij - 1]
                            y2 = olines[id, j, ij]
                            y3 = olines[id, j, ij + 1]
                            tprofs[id, j, i] = pystark.FINTRP(x1, x2, x3, y1, y2, y3, domeg)
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
                            tprof[id, j, i] = pystark.FINTRP(x1, x2, x3, y1, y2, y3, domeg)
                            y1 = olines[id, j, ij_max - 3]
                            y2 = olines[id, j, ij_max - 2]
                            y3 = olines[id, j, ij_max - 1]
                            tprofs[id, j, i] = pystark.FINTRP(x1, x2, x3, y1, y2, y3, domeg)
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
            delta_lambda[i] = wl_centre_angst * delta_omega / (angular_freq_0 + delta_omega)
            # print(delta_lambda[i])
            wprof_nu[i] = (wprof / normal_holtsmark_field) * (2. * np.pi)
            wprofs_nu[i] = (wprofs / normal_holtsmark_field) * (2. * np.pi)
            #        print '%e %e %e %e' %(delta_lambda[i],delta_nu[i],wprof_nu[i],wprofs_nu[i])

        delta_lambda2 = np.concatenate((-delta_lambda[::-1], delta_lambda)) + wl_centre_angst  # + olam0
        delta_nu2 = np.concatenate((-delta_nu[::-1], delta_nu))
        wprof_nu2 = np.concatenate((wprof_nu[::-1], wprof_nu))
        wprofs_nu2 = np.concatenate((wprofs_nu[::-1], wprofs_nu))

        # for some reason, i only get a good agreement with the other models if i take the pure Stark broadened Stehle
        # output and manually convolve it with the Doppler profile -- not sure why...
        ls_sd = wprofs_nu2

        # interpolate onto frequency axis
        ls_sd = np.interp(self.freq_axis, delta_nu2 + self.freq_centre, ls_sd)

        return ls_sd

    def zeeman_split(self, x, x_centre, ls, x_units='Hz'):
        """
         returns input lineshape, with Zeeman splitting accounted for by a simple model
         
        :param x: 
        :param x_centre: 
        :param ls: 
        :param x_units: 
        :return:

        """

        assert x_units in pystark.valid_x_units

        if self.bfield == 0.:
            return ls

        if x_units == 'm':
            freqs, freq_centre, ls = pystark.convert_ls_units(x, x_centre, ls=ls)
        else:
            freqs, freq_centre = x, x_centre

        rel_intensity_pi = np.sin(self.viewangle) ** 2 / 2

        rel_intensity_sigma = 1 / 4 * (1 + np.cos(self.viewangle) ** 2)
        freq_shift_sigma = e / (4 * np.pi * m_e) * self.bfield

        # relative intensities normalised to sum to one

        ls_sigma_minus = rel_intensity_sigma * np.interp(freqs + freq_shift_sigma, freqs, ls)
        ls_sigma_plus = rel_intensity_sigma * np.interp(freqs - freq_shift_sigma, freqs, ls)
        ls_pi = rel_intensity_pi * ls

        return ls_sigma_minus + ls_pi + ls_sigma_plus

    def doppler_convolution(self):
        pass

    def plot(self):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.wl_axis, self.ls_szd)
        ax.axvline(self.wl_centre, color='r')

        # plt.semilogy()
        plt.show()



