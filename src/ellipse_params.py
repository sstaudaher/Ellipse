import copy
import ellipse
import ellipse_io
import ellipse_funcs
import math
import numpy
import pdb
import pickle
import profile_fitting
import matplotlib.pyplot as plt

class Galaxy(object):
    
    """
    
    Calculates the mass of the components of a galaxy based on output
        from the IRAF task ELLIPSE. Also, calculates monte-carlo based
        uncertainties and stores various data about the galaxy.
        
    Required files in ../txt/:
    img_lut.txt, a list of image filenames
    EDGES.dat, a list of general properties for the galaxy from NED
    EDGES_dists.dat, a list of distances from EDD
    EDGES_t.txt, a list of the EDGES morphological type
    EDGES_bars.txt, a list of bar radii
    mark_disk.txt, a list of disk extents 
    
    Required file in ../psf/
    psf_edit.fits, the point-spread-function
    
    Required file in ../img/
    galaxy.fits, the image/mosaic
    
    Required file in ../ellipse/:
    galaxy.txt, text output from IRAF's ELLIPSE
    
    """
    
    def __init__(self, dicts, ellipse_file, psf_file):
        self.set_dicts(dicts)
        self.u = self.__dict__['u']
        self.d25 = self.__dict__['d25']
        self.set_bounds()
        self.n_components = int(len(self.bounds) - 1.0)
        self.set_ellipse(ellipse_file)
        self.set_psf(psf_file)
        self.set_local_params()
        self.ml = 0.5 
        self.ml_unc = 0.1
        self.reach_cal = 280.9
        self.reach_cal_unc = 4.1

    def set_dicts(self, dicts):
        for d in dicts:
            self.__dict__.update(d)

    @property
    def u(self):
        return ellipse_funcs.kpc2u(self.d_kpc)
        
    @u.setter
    def u(self, u):
        self.d_kpc = ellipse_funcs.u2kpc(u)
    
    #D25 is in arseconds    
    @property
    def d25(self):
        return self.r25*2*60.0
    
    #r25 is in arcmin    
    @d25.setter
    def d25(self, d25):
        self.r25 = d25/2.0/60.0
    
    def get_bounds(self):
        bounds = numpy.insert(self.cross_sma, 0, self.r1)
        bounds = numpy.insert(bounds, len(bounds), self.r2)
        return bounds
            
    def set_bounds(self):
        bounds = self.get_bounds()
        self.bounds = bounds
        self.r_break = bounds[1:-1]
                            
    def set_ellipse(self, ellipse_file):
        self.ellipse = ellipse_io.load_ellipse_dict(ellipse_file)

    def set_mosaic(self, mosaic_file):
        self.mosaic, self.mosaic_hdr = ellipse_io.readfits(mosaic_file)

    def set_psf(self, psf_file):
        self.psf, self.psf_hdr = ellipse_io.readfits(psf_file)
        
    def set_img(self):
        self.img, self.hdr = ellipse_io.readfits(self.img_file)
        
    def del_img(self):
        del self.img
        del self.hdr
        
    def set_local_params(self):
        self.arcmin_sma = self.ellipse['SMA']*0.75/60.0
        self.mag, self.mag_lower, self.mag_upper = self.get_mag()
        self.h_local_sma, self.h_local, self.h_local_unc = self.get_h_local()
        self.scale_med = self.get_scale_med()
        self.set_img()
        self.asym = ellipse.asymmetry(self.img, self.ellipse)
        self.del_img()
        
    def get_mag(self):
        mag = ellipse_funcs.intens2mag(self.ellipse['INTENS'])
        mag_lower, mag_upper = ellipse_funcs.mag_unc_solver(self.ellipse['INTENS'], \
                self.ellipse['RMS'], mag)
        return mag, mag_lower, mag_upper
    
    def get_h_local(self):
        h_local_sma, h_local, h_local_unc, h_local_intens, h_local_intens_unc = \
                profile_fitting.local_scale_length_solver( \
                    self.arcmin_sma, self.ellipse['INTENS'], self.ellipse['RMS'])
        return h_local_sma, h_local, h_local_unc
        
    def get_scale_med(self):
        ind = self.ellipse['INTENS'] > 10
        inverse_med = numpy.nanmedian(1.0/self.h_local[ind])
        scale_med = 1.0/inverse_med
        return scale_med
    
    def set_bulge_mass(self):
        self.bulge_mass = self.get_bulge_mass()
            
    def get_bulge_mass(self):
        psf_y_extent, psf_x_extent = numpy.shape(self.psf)
        psf_center = int(psf_y_extent/2)
        psf_normal = self.psf/self.psf[psf_center, psf_center]
        psf_MJy = psf_normal*(0.75/206265.0)**2.0
        psf_bulge = psf_MJy*self.ellipse['INTENS'][0]
        bulge_intens = numpy.sum(psf_bulge)
        return ellipse_funcs.intens2mass(bulge_intens, self.u)
        
    def set_disk_mass(self):
        self.set_disk_params()
        self.set_break_type()
        self.disk_mass = self.get_disk_mass()
     
    def set_disk_params(self):
        self.i0, self.h = self.get_disk_params()
        
    def set_break_type(self):
        self.break_type = self.get_break_type()
        
    def set_disk_params(self):
        self.h, self.h_unc, self.i0, self.i0_unc = \
            profile_fitting.scale_solver(self.ellipse, self.bounds)
    
    def get_break_type(self):
        if len(self.h) == 1:
            break_type = [1]
        else:
            break_type = numpy.zeros(len(self.h) - 1)
            for n in range(len(self.h) - 1):
                if self.h[n] > self.h[n + 1]:
                    if self.r_break[n] < 2.0*self.r_bar: #OLR is <2*r_bar
                        break_type[n] = 2.5  #2.5 is OLR
                    else:
                        break_type[n] = 2    #2 is classical
                if self.h[n] < self.h[n + 1]:
                    break_type[n] = 3
        return break_type
        
    def get_disk_mass(self):
        r_kpc = ellipse_funcs.arcmin2kpc(self.bounds, self.d_kpc)
        h_kpc = ellipse_funcs.arcmin2kpc(self.h, self.d_kpc)
        mass = numpy.zeros(self.n_components)
        for n in range(self.n_components):
            if n == self.n_components - 1:
                if self.break_type[-1] == 3:
                    mass[n] = ellipse_funcs.sersic_mass_int(self.i0[n], h_kpc[n], 200.0, self.d_kpc, 1.0, 1.0)
                else:
                    mass[n] = ellipse_funcs.sersic_mass_partial(self.i0[n], h_kpc[n], self.d_kpc, \
                        1.0, 1.0, r_kpc[n], 200.0)
            else:
                mass[n] = ellipse_funcs.sersic_mass_partial(self.i0[n], h_kpc[n], self.d_kpc, \
                    1.0, 1.0, r_kpc[n], r_kpc[n + 1])
        return mass
    
    def set_monte_mass_unc(self, n_rolls):
        bulge_mass_unc, bulge_mass_frac_unc, disk_mass_unc, disk_mass_frac_unc = \
            self.get_monte_mass_unc(n_rolls)
        self.bulge_mass_unc = bulge_mass_unc
        self.bulge_mass_frac_unc = bulge_mass_frac_unc
        self.disk_mass_unc = disk_mass_unc
        self.disk_mass_frac_unc = disk_mass_frac_unc
        
    def get_monte_mass_unc(self, n_rolls):
        ellipse_dict = copy.deepcopy(self.ellipse)
        u = copy.deepcopy(self.u)
        ml = copy.deepcopy(self.ml)
        reach_cal = copy.deepcopy(self.reach_cal)
        bulge_mass = numpy.zeros(n_rolls)
        disk_mass = numpy.zeros((n_rolls, self.n_components))
        total_mass = numpy.zeros(n_rolls)
        for i in range(n_rolls):
            self.set_rndm_params()
            bulge_mass[i] = self.get_bulge_mass()
            disk_mass[i, :] = self.get_disk_mass()
            total_mass[i] = bulge_mass[i] + sum(disk_mass[i, :])
            self._reset_params(ellipse_dict, u, ml, reach_cal)    
        bulge_mass_unc = numpy.std(bulge_mass)
        bulge_mass_frac_unc = numpy.std(bulge_mass/total_mass)
        disk_mass_unc = numpy.zeros(self.n_components)
        disk_mass_frac_unc = numpy.zeros(self.n_components)
        for i in range(self.n_components):
            disk_mass_unc[i] = numpy.std(disk_mass[:, i])
            disk_mass_frac_unc[i] = numpy.std(disk_mass[:, i]/total_mass)
        self._reset_params(ellipse_dict, u, ml, reach_cal)        
        return bulge_mass_unc, bulge_mass_frac_unc, disk_mass_unc, disk_mass_frac_unc

    def set_rndm_params(self):
        for n in range(len(self.ellipse['INTENS'])):
            self.ellipse['INTENS'][n] = \
                numpy.random.normal(self.ellipse['INTENS'][n], self.ellipse['RMS'][n])
        self.u = numpy.random.normal(self.u, self.u_err)
        self.ml = numpy.random.normal(self.ml, self.ml_unc)
        self.reach_cal = numpy.random.normal(self.reach_cal, self.reach_cal_unc)
        
    def _reset_params(self, ellipse_dict, u, ml, reach_cal):
        self.ellipse = ellipse_dict
        self.u = u
        self.ml = ml
        self.reach_cal = reach_cal
                
def main():
    gals = []
    psf_file = '../psf/psf_edit.fits'
    dicts = ellipse_io.read_dicts()
    for i in range(len(dicts[0]['Galaxy'])):
        my_dicts = ellipse_io.multi_dict_slices(dicts, i)
        ellipse_file = '../ellipse/' + dicts[1]['Galaxy'][i] + '.ch1.txt'
        gal = Galaxy(my_dicts, ellipse_file, psf_file)
        gal.set_bulge_mass()
        gal.set_disk_mass()
        gal.set_monte_mass_unc(1000)
        gals.append(gal)
        print gal.Galaxy
    with open('../pickles/gals.pickle', 'wb') as handle:
      pickle.dump(gals, handle)
    #r_break, hi, ho, break_type = get_break_plot_params(gals)
    #break_stats_plot(r_break, hi, ho, break_type)
    
def get_break_plot_params(gals):
    r_break = []
    hi = []
    ho = []
    break_type = []
    for gal in gals:
        for i in range(len(gal.r_break)):
            r_break.append(gal.r_break[i])
            hi.append(gal.h[i])
            ho.append(gal.h[i+1])
            break_type.append(gal.break_type[i])
        if len(gal.r_break) == 0:
            r_break.append(float('nan'))
            hi.append(float('nan'))
            ho.append(float('nan'))
            break_type.append(gal.break_type[0])
    return r_break, hi, ho, break_type
    
def break_stats_plot(r_break, hi, ho, break_type):
    ellipse_io.init_latex()
    t3 = numpy.asarray(break_type) == 3.0
    t2OLR = numpy.asarray(break_type) == 2.5
    t2 = numpy.asarray(break_type) == 2.0
    t1 = numpy.asarray(break_type) == 1.0
    x = numpy.log10(numpy.asarray(ho)/numpy.asarray(hi))
    y = numpy.asarray(r_break)/numpy.asarray(hi)
    plt.plot(x[t2], y[t2], 'rs', label='Type-2.CT')
    plt.plot(x[t2OLR], y[t2OLR], 'ro', label='Type-2.OLR')
    plt.plot(x[t3], y[t3], 'b^', label='Type-3')
    plt.xlim((-1, 1))
    plt.ylim((0,10))
    plt.xlabel(r'${\rm log}_{10}(h_{\rm o}/h_{\rm i})$')
    plt.ylabel(r'$R_{\rm break}/h_{\rm i}$')
    plt.legend(loc='upper left')
    plt.show()
    pdb.set_trace()
    plt.savefig('break_plot.eps')

if __name__ == '__main__':
    main()