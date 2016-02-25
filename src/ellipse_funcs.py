import numpy
import random
import pdb
from scipy.special import gammainc
from scipy.special import gamma
from math import log10, floor

global VegaZP
global ABZP
VegaZP = 17.6935
ABZP = 20.472

def s_convol(x, y, y_unc, filt):
    x_convol = x[1:-1]
    y_convol = numpy.zeros(len(x_convol))
    y_unc_convol = numpy.zeros(len(x_convol))
    for n in range(len(x_convol)):
        y_convol[n] = sum(y[n: n + 3]*filt)
        y_unc_convol[n] = sum((y_unc[n: n + 3]*filt)**2.0)**0.5
    y_convol = numpy.convolve(y, filt)
    return x_convol, y_convol, y_unc_convol

#Converts magnitudes to mass
def mag2mass(mag):
    MLratio = 0.5
    M_sun = 3.24
    mass = numpy.log10(MLratio) - 0.4*(mag - M_sun - 21.5721)
    return mass

def round_sig(x, sig):
    sig = round(abs(x), sig-int(floor(log10(abs(x))))-1)
    if x < 0:
        sig = -sig
    return sig

def double_profile(I0, h1, h2, R_break, alpha, R):
    exponent = (1.0/alpha*(1.0/h1-1.0/h2))
    S = (1 + numpy.exp(-alpha*R_break))**exponent
    I = S*I0*exp(-R/h1)*(1 + numpy.exp(alpha*(R-R_break)))**exponent
    return I

def gauss_gt_zero(x, x_unc):
    y = random.gauss(x, x_unc)
    while y < 0:
        y = random.gauss(x, x_unc)
    return y

#Calculates the total mass from a section of the surface brightness profile
def sersic_mass_partial(S0, h_kpc, d_kpc, ba, n, r_0_kpc, r_1_kpc):
    mass_inner = sersic_mass_int(S0, h_kpc, r_0_kpc, d_kpc, ba, n)
    mass_outer = sersic_mass_int(S0, h_kpc, r_1_kpc, d_kpc, ba, n)
    mass_total = mass_outer - mass_inner
    return mass_total
    
def Ie2I0(Ie, bn):
    I0 = Ie*numpy.exp(bn)
    return I0
    
def h2Re(h, n):
    Re = h*b**n

def b_calc(n):
    if n == 1:
        b = 1.678
    if n == 4:
        b = 7.669
    return b

def i_solver(i0, h0, h1, r_break):
    m0 = -2.5*numpy.log10(i0)
    m_new = m0 + 2.5/numpy.log(10)*r_break*(1.0/h0 - 1/h1)
    i_new = 10**(-0.4*m_new)
    return i_new

#Calculates a total mass given a sersic profile's parameters
#S0: Surface Brightness in MJy/sr
#h_kpc: scale length in kpc
#d_kpc: distance to galaxy in kpc
#ba: b/a
#n: Sersic index
def sersic_mass_int(S0, h_kpc, R_kpc, d_kpc, ba, n):
    h_rad = h_kpc/d_kpc
    R_rad = R_kpc/d_kpc
    f_int = sersic_int(S0, h_rad, ba, R_rad, n)
    mass_int = intens2mass(f_int, kpc2u(d_kpc))
    return mass_int

def sersic_int(S0, h_rad, ba, R_rad, n):
    f = 2*3.14159*ba*S0*h_rad**2.0*n*gamma(2*n)*gammainc(2*n, (R_rad/h_rad)**1.0/n)
    return f
    
def intens2mass(intens, u):
    M = 0.5*intens*1e6/280.9*10**(0.4*(u + 3.24))
    return M

def gaussian(sigma, a, x):
    f = a*numpy.exp(-x**2.0/2.0/sigma**2.0)
    return f

def kpc2u(d_kpc):
    u = 5*numpy.log10(d_kpc) + 10
    return u

def u2kpc(u):
    kpc = 10**(0.2*(u - 10.0))
    return kpc
    
def u2Mpc(u):
    Mpc = 10**(0.2*u - 5.0)
    return Mpc
    
def kpc_distance(dict, n):
    u_dist = float(dict['u'][n])
    kpc_dist = u2kpc(u_dist)
    return kpc_dist

#Converts to a physical scale [kpc] given
#a distance (also in kpc) from arcmin   
def arcmin2kpc(d_arcmin, dist_kpc):
    d_kpc = 60/206265.0*d_arcmin*dist_kpc
    return d_kpc
    
def arcsec2pc(d_arcsec, dist_kpc):
    d_pc = d_arcsec/206265.0*dist_kpc*1000.0
    return d_pc
    
def intens2mag(intens):
    mag = -2.5*numpy.log10(intens) + VegaZP
    return mag
    
def mag2intens(mag):
    intens = 10**(-0.4*(mag - VegaZP))
    return intens
    
def mag_unc_solver(intens, intens_err, mag):
    if intens_err[0] != intens_err[0]:
        intens_err[0] = intens_err[1]
    intens_upper = intens + intens_err
    intens_lower = intens - intens_err
    ind_lower = numpy.where(intens_lower <= 0.0)
    mag_upper = abs(intens2mag(intens_lower) - mag)
    mag_lower = abs(intens2mag(intens_upper) - mag)
    mag_upper[ind_lower] = 50.0
    return mag_lower, mag_upper