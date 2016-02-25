import numpy
import pyfits
import math
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from ellipse_io import *

"""
A library of functions designed to create residual images from ELLIPSE results

To Do:
Re-write as a class
"""

def asymmetry_solver(gal):
    ellipse_dict = load_ellipse_dict(gal)
    img, hdr = readfits(smooth_dir + gal + '.fits')
    asym = asymmetry(img, ellipse_dict)
    return asym

def lim_calc(ellipse_dict, n):
    sma = ellipse_dict['SMA'][n]
    x0 = ellipse_dict['X0'][n]
    y0 = ellipse_dict['Y0'][n]
    y_lim = [long(round(x0 - sma)), long(round(x0 + sma))]
    x_lim = [long(round(y0 - sma)), long(round(y0 + sma))]
    return x_lim, y_lim

#Finds the total flux within an ellipse
def ellipse_total(stamp, ellipse_dict, n):
    y_extent, x_extent = stamp.shape
    x_center = x_extent/2
    y_center = y_extent/2
    a = ellipse_dict['SMA'][n]
    ellip = ellipse_dict['ELLIP'][n]
    b = a*(1 - ellip)
    theta = ellipse_dict['PA'][n]
    dist_stamp = dist_ellipse_idl(x_extent, y_extent, \
        x_center, y_center, a, b, theta)
    ind = numpy.where(dist_stamp < a)
    total_flux = sum(stamp[ind])
    return total_flux
    
def sma_stamp(stamp, ellipse_dict, n):
    y_extent, x_extent = stamp.shape
    x_center = x_extent/2
    y_center = y_extent/2
    a = ellipse_dict['SMA'][n]
    ellip = ellipse_dict['ELLIP'][n]
    b = a*(1 - ellip)
    theta = ellipse_dict['PA'][n]*3.14159/180.0
    dist_stamp = dist_ellipse_idl(x_extent, y_extent, \
        x_center, y_center, a, b, theta)
    return dist_stamp
    
def ellipse_ind(stamp, ellipse_dict, n):
    ind = numpy.where(dist_stamp > a)
    stamp[ind] = 0
    return stamp
    
def ellipse_stamp(img, ellipse_dict, n):
    x_lim, y_lim = lim_calc(ellipse_dict, n)
    stamp = img[x_lim[0]: x_lim[1], y_lim[0]: y_lim[1]]
    return stamp
    
def ellipse_rot(stamp):
    stamp_180 = numpy.rot90(stamp, 2)
    return stamp_180
    
def asymmetry(img, ellipse_dict):
    asym = [float('nan')]
    for n in range(1, len(ellipse_dict['SMA'])):
        stamp = ellipse_stamp(img, ellipse_dict, n)
        stamp_180 = ellipse_rot(stamp)
        dist_stamp = sma_stamp(stamp, ellipse_dict, n)
        a = ellipse_dict['SMA'][n]
        integrated = stamp*(dist_stamp < a)
        integrated_180 = stamp_180*(dist_stamp < a)
        diff = numpy.nansum(abs(integrated - integrated_180))
        total = numpy.nansum(integrated)
        asym.append(diff/total)
    return asym

#Based on dist_ellipse.pro from W. Landsman.
#Creates an image with pixel values equal to
#the semi-major axis
def dist_ellipse_idl(x_extent, y_extent, \
        x_center, y_center, a, b, theta):
    ratio = float(a)/float(b)
    cosang = math.cos(theta)
    sinang = math.sin(theta)
    x = numpy.arange(x_extent) - x_center
    y = numpy.arange(y_extent) - y_center
    img = numpy.zeros((y_extent, x_extent))
    xcosang = x*cosang
    xsinang = x*sinang
    for n in range(int(y_extent)):
        xfinal = xcosang + y[n]*sinang
        yfinal = -1.0*xsinang + y[n]*cosang
        img[n, :] = ((xfinal*ratio)**2.0 + yfinal**2.0)**0.5
    return img

#I'm not sure I'm doing the grid thing right... but this is
#how you'd do it...
def xy_gen(x_extent, y_extent, x_center, y_center, theta):
    grid = numpy.mgrid[0: y_extent, 0: x_extent]
    y = grid[0] 
    x = grid[1]
    x_rot, y_rot = rotation(x, y, theta)
    return x_rot, y_rot
    
def ell_index(x_rot, y_rot, a, b):
    ellip_value = ellipse_eq(x_rot, y_rot, a, b)
    x_ind, y_ind = numpy.where(ellip_value <= 1)
    return x_ind, y_ind
    
#Calculates the index of x and y within rad pixels
def ellipse_index_radius(x_extent, y_extent, \
        x_center, y_center, a, b, theta):
    x_ind = []
    y_ind = []
    for i in range(x_extent):
        for j in range(y_extent):
            x, y = rotation(i, j, theta)
            ellip_value = ellipse_eq(x, y, a, b)
            if ellip_value <= 1:
                x_ind.append(x)
                y_ind.append(y)
    return [x_ind, y_ind]

def ellipse_uncertainty(img, x_extent, y_extent):
    for n in range(a):
        outer_rad = a[n]+0.5
        inner_rad = a[n]-0.5
        [x_ind_outer, y_ind_outer] = \
            ellipse_index_radius(x_extent, y_extent, outer_rad, n)
        [x_ind_inner, y_ind_inner] = \
            ellipse_index_radius(x_extent, y_extent, inner_rad, n)
        outer = img[y_ind_outer, x_ind_outer]
        inner = img[y_ind_outer, x_ind_outer]

#Rotate an x and y matrix
def rotation(x, y, theta):
    x_rot = x*math.cos(theta) - y*math.sin(theta)
    y_rot = x*math.sin(theta) + y*math.cos(theta)
    return x_rot, y_rot

#Finds the value of the equation for an ellipse
def ellipse_eq(x, y, a, b):
    ellipse_value = x**2.0/a**2.0 + y**2.0/b**2.0
    return ellipse_value

#From Ramanujan, accurate to 5th order
def ellipse_circumference(a, b):
    x = (a - b)/(a + b)
    pi = 3.14159
    C = pi*(a + b)*(1 + 3*x**2/(10 + (4 - 3*x**2)))
    return C

#The parametric form, which generates x and y
def ellipse_parametric(a, b, theta):
    x = a*math.cos(theta)
    y = b*math.sin(theta)
    return x, y

#Creates an interpolated image given a number of ellipses        
def make_image(x_extent, y_extent, \
    intens, intens_err, x_center, y_center, a, b, theta):
    img = numpy.zeros((y_extent, x_extent))
    for n in range(len(intens)):
        n_circles = int(2.0*3.14159*a[n])
        for m in range(n_circles):
            angle = m*2*3.14159/n_circles
            x, y = ellipse_parametric(a[n], b[n], angle)
            x_rot, y_rot = rotation(x, y, (theta[n]+90)*3.14159/180.0)
            x_centered = x_rot + x_center[n]
            y_centered = y_rot + y_center[n]
            x_rounded = round(x_centered)
            y_rounded = round(y_centered)
            if (x_rounded >= 0) and (x_rounded < x_extent) and \
                (y_rounded >= 0) and (y_rounded < y_extent):
                    img[y_rounded, x_rounded] = intens[n]
    return img
    
#Solves a linear equation given two 2-d points
def linear_solver(x1, x2, y1, y2):
    m = (y1 - y2)/(x1 - x2)
    b = y1 - m*x1
    return m, b
        
#Determines where a transition (from a number to 0 occurs)
#Specifically, this is calibrated for interpolating values
#betweeen two ellipses
def transition(array):
    transition = []
    #I want to start it when we go from a value to a zero
    starter = 0
    if array[0] == 0:
        switch = 0
    else:
        switch = 1
    for n in range(1, len(array)):
        if array[n] == 0 and switch == 1:
            starter = 1
            switch = 0
            transition.append(n-1)
        if array[n] != 0 and switch == 0:
            switch = 1
            if starter == 1:
                transition.append(n)
    #If it's odd that means it's something like 11000 at end, so
    #we must remove the final transition
    if len(transition) % 2 == 1:
        transition.pop(-1)
    return transition

#Calculates the nearest distance to x and y
def nearest_dist(x, y, x_arr, y_arr):
    dist = (x - x_arr)**2 + (y - y_arr)**2
    print min(dist)
    x_ind, y_ind = numpy.where(dist == min(dist))
    x_ind = x_ind[0]
    y_ind = y_ind[0]
    return x_arr[x_ind], y_arr[y_ind]
    
def nearest_neighbors(a_img, a, ellip_img):
    final_img = numpy.array(ellip_img)
    x_ind, y_ind = numpy.where(ellip_img != 0)
    logic_img = numpy.logical_and(a_img < a[-1], ellip_img == 0)
    x0_ind, y0_ind = numpy.where(logic_img)
    for n in range(len(x0_ind)):
        x = x0_ind[n]
        y = y0_ind[n]
        x_nearest, y_nearest = nearest_dist(x, y, x_ind, x_ind)
        intens = ellip_img[x_nearest, y_nearest]
        logic_img = numpy.logical_and(ellip_img != intens, ellip_img != 0)
        x_second_ind, y_second_ind = numpy.where(logic_img)
        x_second_nearest, y_second_nearest = nearest_dist(x, y, \
            x_second_ind, y_second_ind)
        x1 = 0.0
        x2 = ((x_nearest - x_second_nearest)**2 + \
            (y_nearest - y_second_nearest)**2)**0.5
        x3 = sqrt((x_nearest - x)**2 + (y_nearest - y)**2)
        intens_second = ellip_img[x_second_nearest, y_second_nearest]
        m, b = linear_solver(x1, x2, math.log(intens), math.log(intens_second))
        final_intens = math.exp(m*x3 + b)
        final_img[x, y] = final_intens
    return final_img
    
def fill_in_array(array, transition):
    for n in range(0, len(transition), 2):
        x1 = transition[n]
        x2 = transition[n+1]
        y1 = array[x1]
        y2 = array[x2]
        if y1 != y2:
            y1 = math.log(array[x1])
            y2 = math.log(array[x2])
            m, b = linear_solver(x1, x2, y1, y2)
            for i in range(x1 + 1, x2):             
                array[i] = math.exp(m*i + b)      
    return array
        
def bound_calculator(extent, center, dist_max):
    if center - dist_max > 0:
        lower_bound = center - dist_max
    else:
        lower_bound = 0
    if center + dist_max < extent - 1:
        upper_bound = center + dist_max
    else:
        upper_bound = extent
    return lower_bound, upper_bound
        
def fill_in_image(img, x_extent, y_extent, \
        intens, intens_err, x_center, y_center, a, b, theta):
    final_img = numpy.zeros((y_extent, x_extent))
    dist_max = int(a[-1])
    x0 = int(x_center[-1])
    y0 = int(y_center[-1])
    x_lower_bound, x_upper_bound = \
        bound_calculator(x_extent, x0, dist_max)
    y_lower_bound, y_upper_bound = \
        bound_calculator(y_extent, y0, dist_max)
    #Loops through possible non-zero values
    for i in range(x_lower_bound, x_upper_bound):
        column = img[:, i]
        trans = transition(column)
        filled_array = fill_in_array(column, trans)
        final_img[:, i] = filled_array
    for j in range(y_lower_bound, y_upper_bound):
        row = final_img[j, :]
        trans = transition(row)
        filled_array = fill_in_array(row, trans)
        final_img[j, :] = filled_array
    return final_img
    
#Finds the final values of the ellipse paramaters and the difference
#between the final two a and b values
def end_values(x_extent, y_extent, \
        x_center, y_center, a, b, theta):
    x_cen = x_center[-1]
    y_cen = y_center[-1]
    a_initial = a[-1]   #We start where ellipse left off
    b_initial = b[-1]
    angle = theta[-1]*3.14159/180
    a_diff = a[-1] - a[-2]
    b_diff = b[-1] - b[-2]  
    return x_cen, y_cen, a_initial, b_initial, angle, a_diff, b_diff  

def closest_side(x_extent, y_extent, a_img):
    left = numpy.min(a_img[:, 0])
    right = numpy.min(a_img[:, x_extent - 1])
    up = numpy.min(a_img[y_extent - 1, :])
    down = numpy.min(a_img[0, :])
    return min(left, right, up, down)
   
#Continues calculating values to the end of the image
def ellipse_continue(img, x_extent, y_extent, \
        intens, intens_err, x_center, y_center, a, b, theta, rad):
    x_cen, y_cen, a_initial, b_initial, angle, a_diff, b_diff = \
        end_values(x_extent, y_extent, x_center, y_center, a, b, theta)
    a_img = dist_ellipse_idl(x_extent, y_extent, \
        x_cen, y_cen, a_initial, b_initial, angle)
    a_remaining = closest_side(x_extent, y_extent, a_img) - a_initial
    num_measurements = long(a_remaining/a_diff) - 1
    if num_measurements > 0:
        for n in range(num_measurements):
            a_current = a_initial + n*a_diff
            logic_img = numpy.logical_and(a_img < a_current + rad, \
                a_img > a_current - rad)
            y_ind, x_ind = numpy.where(logic_img)
            annulus = img[y_ind, x_ind]
            nan_ind = numpy.where(annulus == annulus)
            annulus = annulus[nan_ind]
            clipped_annulus, lower, upper = sigmaclip(annulus)
            if n == 0:
               intens_diff = numpy.mean(clipped_annulus) - intens[-1]
            else:
               intens = numpy.append(intens, numpy.mean(clipped_annulus) - intens_diff)
               intens_err = numpy.append(intens_err, numpy.std(clipped_annulus))
               x_center = numpy.append(x_center, x_cen)
               y_center = numpy.append(y_center, y_cen)
               a = numpy.append(a, a_current)
               b = numpy.append(b, a_current*b_initial/a_initial)
               theta = numpy.append(theta, angle*180/3.14159)
    return intens, intens_err, x_center, y_center, a, b, theta
         
#The main program, creates an image from the ellipse results
def ellipse(x_extent, y_extent, \
    intens, intens_err, x_center, y_center, a, b, theta):
    img = make_image(img_filename, x_extent, y_extent, \
        intens, intens_err, x_center, y_center, a, b, theta)
    filled_img = fill_in_image(img, x_extent, y_extent, \
        intens, intens_err, x_center, y_center, a, b, theta)
    return filled_img
    
#Given a img and ellipse file, calculates a residual image
def ellipse_resid(img_file, ellip_file):
    img = readfits(img_file)
    y_extent, x_extent = img.shape
    rad = 2.0
    intens, intens_err, x_center, y_center, a, b, theta = \
        read_ellipse(ellip_file)
    intens, intens_err, x_center, y_center, a, b, theta = \
        ellipse_continue(img, x_extent, y_extent, \
            intens, intens_err, x_center, y_center, a, b, theta, rad)
    ellip_img = make_image(x_extent, y_extent, \
        intens, intens_err, x_center, y_center, a, b, theta)
    #writefits('test.fits',ellip_img)
    filled_img = fill_in_image(ellip_img, x_extent, y_extent, \
        intens, intens_err, x_center, y_center, a, b, theta)
    #writefits('filled_test.fits', filled_img)
    resid_img = img - filled_img
    return resid_img

#Calculates then writes a residual image from ellipse output
def resid(img_file, ellip_file):
    resid_img = ellipse_resid(img_file, ellip_file)
    writefits('resid.fits', resid_img)

def main():
    cat_file = '../txt/mosaic_lut.txt'
    mosaic_lut = loadfile(cat_file)
    for n in range(len(mosaic_lut)):
        gal, mosaic = mosaic_lut[n].split()
        print gal
        img_file = '../smoothed/'+mosaic
        ellip_file = '../ellipse_results/'+gal+'.ch1.txt'
        resid_img = ellipse_resid(img_file, ellip_file)
        writefits('../resids/'+gal+'.fits', resid_img)

if __name__ == '__main__':
    main()