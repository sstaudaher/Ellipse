import numpy
import math
from scipy.optimize import curve_fit
from lmfit import Model
import pdb
import matplotlib.pyplot as plt
import ellipse_funcs

def scale_length_fit_nogal(arcmin_sma, intens, intens_err):
    if len(intens) < 2:
        popt = [float('nan'), float('nan')]
        perr = [float('nan'), float('nan')]
    else:
        try:
            popt, pcov = curve_fit(scale_length_func, arcmin_sma, intens, \
                p0 = [1.0, 2.0], sigma = intens_err, absolute_sigma = False)
            perr = numpy.sqrt(numpy.diag(pcov))
        except RuntimeError:
            popt = [float('nan'), float('nan')]
            perr = [float('nan'), float('nan')]   
    return popt, perr
    
def scale_solver(ellipse_dict, bounds):
    n_points = len(bounds) - 1
    h = numpy.zeros(n_points)
    h_unc = numpy.zeros(n_points)
    i0 = numpy.zeros(n_points)
    i0_unc = numpy.zeros(n_points)
    for n in range(len(bounds) - 1):
        arcmin_sma = ellipse_dict['SMA']*0.75/60.0
        ind = numpy.where((arcmin_sma > \
            bounds[n])*(arcmin_sma < bounds[n + 1]))
        fit, err = scale_length_fit_nogal(arcmin_sma[ind[0]], \
            ellipse_dict['INTENS'][ind[0]], ellipse_dict['RMS'][ind[0]])
        h[n], h_unc[n], i0[n], i0_unc[n] = fit[0], err[0], fit[1], err[0]
    return h, h_unc, i0, i0_unc

#a[0] < a[1]
#b[0] < b[1]
def overlap_test(a, b):
    flag = 1
    if a[0] > b[0] and a[0] > b[1]:
        flag = 0
    if a[1] < b[0] and a[1] < b[1]:
        flag = 0
    return flag
    
def overlap_calc(h, h_unc):
    h_unc = 0.5*h_unc
    overlap = numpy.zeros(len(h) - 1)
    for n in range(len(h) - 1):
        range_1 = ((h[n] - h_unc[n]), (h[n] + h_unc[n]))
        range_2 = ((h[n + 1] - h_unc[n + 1]), (h[n + 1] + h_unc[n + 1]))
        overlap[n] = overlap_test(range_1, range_2)
    return overlap

#r_break algorithm:
#1.) Define the boundaries
#2.) Calculate the scale length for each section
#3.) If the scale lengths don't differ by X-sigma, remove a break
#4.) Recalculate again
#5.) See if they are different by X-sigma again
def break_solver(ellipse_dict, bounds):
    arcmin_sma = ellipse_dict['SMA']*0.75/60.0
    h = numpy.zeros(len(bounds) - 1)
    h_unc = numpy.zeros(len(bounds) - 1)
    for n in range(len(bounds) - 1):
        ind = numpy.where((arcmin_sma > \
            bounds[n])*(arcmin_sma < bounds[n + 1]))
        if len(ind[0]) < 2:
            fit, err = [float('inf'), float('inf')], [float('inf'), float('inf')]
        else:
            fit, err = scale_length_fit_nogal(arcmin_sma[ind[0]], \
                ellipse_dict['INTENS'][ind[0]], ellipse_dict['RMS'][ind[0]])
        h[n] = fit[0]
        h_unc[n] = err[0]
    return h, h_unc
    
def bounds_inf(ellipse_dict, bounds):
    flag = 1
    while flag == 1:
        h, h_unc = break_solver(ellipse_dict, bounds)
        ind = numpy.where(h_unc == float('inf'))
        if len(ind[0]) == 0:
            flag = 0
        else:
            bounds = numpy.delete(bounds, ind[0][0] + 1)
    return bounds, h, h_unc
    
def breaks(ellipse_dict, bounds):
    flag = 1
    while flag == 1:
        bounds, h, h_unc = bounds_inf(ellipse_dict, bounds)
        overlap = overlap_calc(h, h_unc)
        ind = numpy.where(overlap == 1)
        if len(ind[0]) == 0:
            flag = 0
        else:
            bounds = numpy.delete(bounds, ind[0][0] + 1)
            if len(bounds) == 2:
                flag = 0
    return bounds
        
def disk_disk_func(R, I0, h1, h2, R_break):
    I = numpy.zeros(len(R))
    I1 = I1_calc(I0, h1, h2, R_break)
    for n in range(len(R)):
        if R[n] < R_break:
            I[n] = sersic_intens_I0(R[n], I0, h1, 1.0)
        else:
            I[n] = sersic_intens_I0(R[n], I1, h2, 1.0)
    return I

def multi_disk_fit(arcmin_sma, intens, intens_err, bounds, psf):
    h_mat_0 = numpy.ones(len(bounds - 1))
    mod = Model(multi_disk_func)
    pars = mod.make_params(I0 = 1.0, h_mat = h_mat_0, bounds = bounds, PSF = PSF)
    pars['bounds'].vary = False
    pars['PSF'].vary = False
    results = mod.fit(intens, params = pars, x = arcmin_sma, weights = 1.0/intens_err)
    pdb.set_trace()

def multi_disk_func(R, I0, h_mat, bounds, psf):
    bounds = [0, 3.0, 20.0]
    n_points = long(bounds[-1]*60/0.75/4.0)
    
    model_x = numpy.arange(n_points)*4.0*0.75/60.0
    S = numpy.zeros(n_points)
    for n in range(len(h_mat)):
        logic = (bounds[n] <= model_x)*(model_x <= bounds[n + 1])
        ind = numpy.where(logic)
        if n == 0:
            I = I0
        else:
            I = ellipse_funcs.i_solver(I, h_mat[n - 1], h_mat[n], bounds[n])
        S[ind] = sersic_intens_disk(model_x[ind], I, h_mat[n])
    S[0] = 100.0
    S_convol = numpy.convolve(S, psf)
    S_convol = S_convol[0:n_points]
    pdiff = abs(1 - S/S_convol)
    plt.plot(model_x, 2.5*numpy.log10(S_convol))
    plt.plot(model_x, 2.5*numpy.log10(psf/psf[0]*S_convol[0]))
    plt.show()
    pdb.set_trace()
    return S

def multi_fit(R, I, bounds, psf):
    popt, pcov = curve_fit(lambda bounds, psf: multi_disk_func(R, I0, h_mat, bounds, psf), \
        R, I, p0 = p0, sigma = I_unc)
    popt, pcov = curve_fit(disk_disk_func, arcmin_sma, \
        intens, p0 = p0, sigma = intens_err)

def disk_halo_func(R, I0, h, S0):
    disk = sersic_intens_I0(R, I0, h, 1.0)
    halo = halo_u_model(R, S0)
    return disk + halo

def halo_u_model(R_kpc, S0):
    alpha = 1.26
    a_h = 5.2
    R0 = 30.0
    S = S0*((1.0 + (R0/a_h)**2.0)/(1.0 + (R_kpc/a_h)**2.0))**alpha
    return S

def disk_disk_fit(arcmin_sma, intens, intens_err):
    p0 = [1.0, 2.0, 3.0, 0.5*max(arcmin_sma)]
    popt, pcov = curve_fit(disk_disk_func, arcmin_sma, \
        intens, p0 = p0, sigma = intens_err)
    perr = numpy.sqrt(numpy.diag(pcov))
    return popt, perr
    
def disk_halo_fit(kpc_sma, intens, intens_err):
    p0 = [1.0, 1.0, 0.0015]
    popt, pcov = curve_fit(disk_halo_func, kpc_sma, \
        intens, p0 = p0, sigma = intens_err)
    perr = numpy.sqrt(numpy.diag(pcov))
    return popt, perr
    
def fit_solver(ellipse_dict, kpc_dist):
    arcmin_sma = ellipse_dict['SMA']*0.75/60.0
    kpc_sma = arcmin_sma/60.0/206265.0*kpc_dist
    disk_disk_params, disk_disk_err = disk_disk_fit(arcmin_sma, \
        ellipse_dict['INTENS'], ellipse_dict['RMS'])
    disk_halo_params, disk_halo_err = disk_halo_fit(kpc_sma, \
        ellipse_dict['INTENS'], ellipse_dict['RMS'])
    return [disk_disk_params, disk_disk_err], [disk_halo_params, disk_halo_err]
        
def two_disk_profile(m1, m2, x_break, b1, R):
    S = numpy.zeros(len(R))
    b2 = b2_calc(m1, m2, x_break, b1)
    for n in range(len(R)):
        if R < x_break:
            S[n] = m1*R + b1
        else:
            S[n] = m2*R + b2
    return S

def I1_calc(I0, h1, h2, R_break):
    m1 = 2.5/math.log(10)/h1
    m2 = 2.5/math.log(10)/h2
    b1 = -2.5*math.log10(I0)
    b2 = b2_calc(m1, m2, R_break, b1)
    I1 = 10.0**(-0.4*b2)
    return I1

def b2_calc(m1, m2, x_break, b1):
    b2 = (m1 - m2)*x_break + b1
    return b2

def linear_func(x, a, b):
    return a*x + b
        
def b_solver(n):
    if n == 4:
        b = 7.669
    if n == 1:
        b = 1.678
    #From Capaccioli89
    if n != 4 and n != 1 and \
        0.5 < n and n < 10:
        b = 1.9992*n - 0.3271
    if n >= 10:
        b = 2.0*n - 1.0/3.0
    return b
        
def Ie_to_I0(Ie, n):
    b = b_solver(n)
    I0 = Ie*numpy.exp(b)
    return I0
    
def h_to_Re(h):
    if n == 1:
        Re = 1.678*h
    if n == 4:
        Re = 3466*h
    return Re
    
def sersic_intens_I0(R, I0, h, n):
    I = I0*numpy.exp(-1.0*(R/h)**(1.0/n))
    return I
    
def sersic_intens_disk(R, I0, h):
    I = I0*numpy.exp(-1.0*(R/h))
    return I
    
def sersic_mag_u0(R, u0, h, n):
    mag = u0 + 2.5/numpy.log(10)*(R/h)**(1.0/n)
    return mag
    
def sersic_intens(R, Ie, Re, n):
    b = b_solver(n)
    I = Ie*numpy.exp(-b*(R/Re)**(1.0/n) - 1.0)
    return I

def sersic_mag(R, ue, Re, n):
    b = b_solver(n)
    u = ue + 2.5*b/numpy.log(10)*((R/Re)**(1.0/n) - 1.0)
    return u
    
def scale_length_func(x, a, b):
    return b*numpy.exp(-x/a)

#Tests to find if a rise is too big for curve_fit to handle
def rise_calc(y):
    rise = y[0:-1] - y[1:]
    ind = numpy.where(rise < 0)
    if len(ind[0]) > 0:
        med_rise = numpy.median(rise)
        ind_med = numpy.where(abs(rise[ind]) > 2*med_rise)
        if len(ind_med[0]) > 0:
            flag = 0
        else:
            flag = 1
    else:
        flag = 1
    return flag

def scale_length_fit(gal, arcmin_sma, intens, intens_err):
    try:
        popt, pcov = curve_fit(scale_length_func, arcmin_sma, intens, \
            p0 = [1.0, 2.0], sigma = intens_err, absolute_sigma = False)
        perr = numpy.sqrt(numpy.diag(pcov))
    except RuntimeError:
        popt = [float('nan'), float('nan')]
        perr = [float('nan'), float('nan')]
    return popt, perr
    
def local_scale_length_solver(gal, arcmin_sma, intens, intens_err):
    n_points = 11
    n_array = len(intens) - n_points + 1
    if n_array > n_points:
        ind = numpy.arange(n_points)
        scale_length_sma = numpy.zeros(n_array)
        scale_length = numpy.zeros(n_array)
        scale_length_unc = numpy.zeros(n_array)
        scale_intens = numpy.zeros(n_array)
        scale_intens_unc = numpy.zeros(n_array)
        for n in range(n_array):
            popt, perr = scale_length_fit(gal, arcmin_sma[ind], intens[ind], intens_err[ind])
            scale_length_sma[n] = numpy.mean(arcmin_sma[ind])
            scale_length[n] = popt[0]
            scale_length_unc[n] = perr[0]
            scale_intens[n] = popt[1]
            scale_intens_unc[n] = perr[1]
            ind = ind + 1
        #ind = numpy.where((scale_length_sma > 0.25)* \
        #                  (abs(scale_length_unc/scale_length) < 100.0))
        #scale_length_sma = scale_length_sma[ind]
        #scale_length = scale_length[ind]
        #scale_length_unc = scale_length_unc[ind]
        #scale_intens = scale_intens[ind]
        #scale_intens_unc = scale_intens_unc[ind]
        return scale_length_sma, scale_length, scale_length_unc, \
            scale_intens, scale_intens_unc
    else:
        pdb.set_trace()
        print "Error: n_points < n_array"