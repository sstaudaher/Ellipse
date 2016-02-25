import numpy
import pyfits
import math
import matplotlib
import matplotlib.pyplot as plt

#A library of input/output functions designed to work with IRAF
#ELLIPSE outputs

def multi_dict_slices(my_dicts, n):
    dict_slices = []
    for my_dict in my_dicts:
        dict_slices.append(dict_slice(my_dict, n))
    return dict_slices
        
#Returns a single row in a dictionary
def dict_slice(my_dict, n):
    d = {}
    keys = my_dict.keys()
    for key in my_dict.keys():
        d[key] = my_dict[key][n]
    return d

#Reads in a bunch of text files into dictionaries    
def read_dicts():
    #gal, mosaic_file = mosaic_dictionary['Galaxy'], mosaic_dictionary['Mosaic']
    dirs = set_dirs()
    edges_dictionary = ellipse_io.float_readcol(dirs['txt_dir'] + 'EDGES.dat', 0, 2)
    edges_distances = ellipse_io.float_readcol(dirs['txt_dir'] + 'EDGES_dists.txt', 1, 2)
    disk_dictionary = ellipse_io.float_readcol(dirs['txt_dir'] + 'mark_disk.txt', 0, 1)
    edges_bars = ellipse_io.float_readcol(dirs['txt_dir'] + 'EDGES_bars.txt', 0, 1)
    edges_t = ellipse_io.float_readcol(dirs['txt_dir'] + 'EDGES_t.txt', 1, 2)
    cross_dict = readcross()
    return [edges_dictionary, edges_distances, disk_dictionary, edges_bars, edges_t, cross_dict]

def set_dirs():
    dirs = {}
    dirs['txt_dir'] = '../txt/'
    dirs['mosaic_dir'] = '../imgs/'
    dirs['region_dir'] = '../regions/'
    dirs['sex_results_dir'] = '../se_catalogs/'
    dirs['smooth_dir'] = '../smoothed/'
    dirs['ellipse_dir'] = '../ellipse/'
    dirs['psf_dir'] = '../psf/'
    return dirs

def readcross():
    lines = ellipse_io.loadfile('../txt/s_breaks.txt')
    cross_dict = {}
    gal = []
    cross_sma = []
    for i in range(len(lines) - 1):
        tmp = lines[i].split()
        gal.append(tmp[0])
        cross_sma.append([float(i) for i in tmp[1:]])
    cross_dict['Galaxy'] = gal
    cross_dict['cross_sma'] = cross_sma
    return cross_dict

def ellipse_region(region_file, ellipse_dict):
    a, b = get_a_b(ellipse_dict)
    ellipse_ds9(region_file, ellipse_dict['X0'], ellipse_dict['Y0'], \
        a, b, ellipse_dict['PA'])

def r25_last_ellipse(filename, ellipse_dict, edges_dictionary, n):
    a, b = get_a_b(ellipse_dict)
    edges_dictionary['d25']
    
def reg_23(filename, ellipse_dict):
    mag_23 = abs(ellipse_dict['mag'] - 23)
    ind = numpy.where(min(mag_23) == mag_23)
    a, b = get_a_b(ellipse_dict)
    ellipse_ds9(filename, ellipse_dict['X0'][ind], ellipse_dict['Y0'][ind], \
                a[ind], b[ind], ellipse_dict['PA'][ind])

#Creates a website to show off your sb_plots
def web_maker(gals):
    with open('profiles.html', 'w') as f:
        f.write('<html>\n')
        f.write('<title>Surface Brightness Profiles</title>\n')
        f.write('<body>\n')
        f.write('\n')
        for gal in gals:
            f.write('<tr>\n')
            f.write('  <td>' + gal +'<br><img src="images/' + gal + \
                    '.jpeg" heights="400" width="500"></td>\n')
            f.write('  <td><img src="scale_length/' + gal + \
                    '.jpeg" height="400" width = "500"></td><br>\n')
            #f.write('  <td><img src="fixed/scale_length/' + gal + \
            #        '.jpeg" height="400" width = "500"></td><br>\n')
            f.write('</tr>\n')
            f.write('\n')

def ds9_test(img_names, region_names):
    file = open('ds9_regions.csh', 'w')
    for n in range(len(img_names)):
        file.write('ds9 ' + img_names[n] + '-regions ' + 
            regions_names[n] + '\n')
    file.close()

#Sets latex up
def init_latex():
    plt.rc('text', usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amssymb}']

def concat(mat):
    new = []
    for n in range(len(mat)):
        if type(mat[n]) is numpy.float64:
            new.append(mat[n])
        else:
            for i in range(len(mat[n])):
                new.append(mat[n][i])
    return numpy.array(new)

def two_concat(a, b):
    new = numpy.zeros(len(a) + len(b))
    new[0:len(a)] = a
    new[len(a): len(a) + len(b)] = b
    return new

def float_convert(mylist):
    if type(mylist) is list:
        mylist = numpy.asarray(mylist).astype(float)
    else:
        print 'This is not a list =/'
    return mylist

#"Untangles" nested lists
def list_untangle(mylist):
    untangled = []
    for i in range(len(mylist)):
        if type(mylist[i]) is float:
            untangled.append(mylist[i])
        else:
            for j in range(len(mylist[i])):
                untangled.append(mylist[i][j])
    return untangled

#Adds multiple variables to a dictionary
def ellipse_dict_add_var(ellipse_dict, var_name, var):
    for n in range(len(var_name)):
        ellipse_dict[var_name[n]] = var[n]
    return ellipse_dict

#Updates the ellipse dictionary with a new index
def ellipse_dict_ind_update(ellipse_dict, ind):
    for var in ellipse_dict:
        ellipse_dict[var] = ellipse_dict[var][ind]
    return ellipse_dict

def load_ellipse_dict(filename):
    ellipse_dict = readcol(filename, 2, 5)
    ellipse_dict = ellipse_dict_fix(ellipse_dict)
    return ellipse_dict

#Changes INDEF to NaN and converts
#all strings to floats
def ellipse_dict_fix(ellipse_dict):
    for hdr in ellipse_dict:
        for n in range(len(ellipse_dict[hdr])):
            ellipse_dict[hdr][n] = indef_fix(ellipse_dict[hdr][n])
            ellipse_dict[hdr][n] = float(ellipse_dict[hdr][n])
        ellipse_dict[hdr] = numpy.zeros(len(ellipse_dict[hdr])) + \
            ellipse_dict[hdr]
    ellipse_dict['X0'] = ellipse_dict['X0'] - 1
    ellipse_dict['Y0'] = ellipse_dict['Y0'] - 1
    if ellipse_dict['RMS'][0] != ellipse_dict['RMS'][0]:
        ellipse_dict['RMS'][0] = ellipse_dict['RMS'][1]
    return ellipse_dict

def hash_fix(hdr):
    if hdr[0] == '#':
        hdr = hdr[1:]
    return hdr

#returns u and u_unc given the distance dictionary
def u_finder(dictionary, n):
    u = float(dictionary['u'][n])
    u_unc = float(dictionary['u_err'][n])
    return u, u_unc

#Filename: a string containing the file's name
#hdrline: the line that the header appears (starting at 0)
#startline: the line that the data starts (starting at 0)
def readcol(filename, hdrline, startline):
    dictionary = {}
    lines = loadfile(filename)
    hdr = lines[hdrline].split()
    hdr = hash_fix(hdr)
    for i in range(len(hdr)):
        dictionary[hdr[i]] = []
    for n in range(startline, len(lines)):
        column = lines[n].split()
        if column != []:
            for i in range(len(hdr)):
                dictionary[hdr[i]].append(column[i])
    return dictionary
    
def float_readcol(filename, hdrline, startline):
    dictionary = readcol(filename, hdrline, startline)
    float_dictionary = float_dict(dictionary)
    return float_dictionary
    
def float_dict(dictionary):
    new_dict = {}
    for key in dictionary.keys():
        try:
            new_dict[key] = [float(i) for i in dictionary[key]]
        except ValueError:
            new_dict[key] = dictionary[key]
    return new_dict
    
def string_list_space(list):
    string = ''
    for n in range(len(list)):
        string = string + str(list[n]) + ' '
    return string
    
def writedict(dictionary, filename):
    dict_keys = dictionary.keys()
    with open(filename, 'w') as f:
        string = ''
        for n in range(len(dict_keys)):
            string = string + dict_keys[n] + ' '
        string = string + '\n'
        f.write(string)
        for n in range(len(dictionary[dict_keys[0]])):
            string = ''
            for i in range(len(dict_keys)):
                string = string + str(dictionary[dict_keys[i]][n]) + ' '
            string = string + '\n'
            f.write(string)

def mosaic_lut_extract(lines):
    gal = []
    mosaic_file = []
    for n in range(1, len(lines)):
        if lines[n] != []:
            mosaic_lut_split = lines[n].split()
            gal.append(mosaic_lut_split[0])
            mosaic_file.append(mosaic_lut_split[1])
    return gal, mosaic_file

#Everything in pixels or degrees
def ellipse_ds9(filename, x_center, y_center, a, b, theta):
    with open(filename, 'w') as f:
        for n in range(len(x_center)):
            string = 'ellipse(' + str(x_center[n] + 1) + ',' + str(y_center[n] + 1) + ',' + \
                     str(a[n]) + ',' + str(b[n]) + ',' + str(theta[n] + 90) + ')\n'
            f.write(string)
                        
#Reads a fit file, simply...
def readfits(filename):
    with pyfits.open(filename) as hdulist:
        scidata = hdulist[0].data
        hdr = hdulist[0].header
    return scidata, hdr

#Writes a fit file, simply...
def writefits_nohdr(filename, img):
    hdu = pyfits.PrimaryHDU(img)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)
    hdulist.close()

#Writes a fit file, simply...
def writefits(filename, img, hdr):
    hdu = pyfits.PrimaryHDU(img, header = hdr)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)
    hdulist.close()

#loads the file into a list split by line breaks
def loadfile(filename):
    file = open(filename, 'r')
    everyline = file.read()
    lines = everyline.split('\n')
    file.close()
    return lines

#Finds what row the data begins on by
#skipping empty rows and commented rows
def header_row(lines):
    for n in range(0, len(lines)):
        split_line = lines[n].split()
        empty_line_test = split_line != []
        if empty_line_test:
            comment_test = split_line[0] != '#'
            if comment_test:
                header_row = n
                break
    return header_row

#Finds and returns the column number for a specific variable
def find_col(lines, colname):
    for n in range(0, len(lines)):
        split = lines[n].split()
        for col in range(len(split)):
            if split[col] == colname:
                return col - 1
#Sets latex up
def init_latex():
    plt.rc('text', usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amssymb}']

def init_ellipse_params(num_lines):
    intens = numpy.zeros(num_lines)
    intens_err = numpy.zeros(num_lines)
    x_center = numpy.zeros(num_lines)
    y_center = numpy.zeros(num_lines)
    a = numpy.zeros(num_lines)
    b = numpy.zeros(num_lines)
    theta = numpy.zeros(num_lines)
    return intens, intens_err, x_center, y_center, a, b, theta

#Determines what column the ellipse paramters are located in  
def ellipse_param_inds(lines):
    intens_ind = find_col(lines, "INTENS")
    intens_err_ind = find_col(lines, "RMS")
    x_center_ind = find_col(lines, "X0")
    y_center_ind = find_col(lines, "Y0")
    a_ind = find_col(lines, "SMA")
    ellipticity_ind = find_col(lines, "ELLIP")
    theta_ind = find_col(lines, "PA")
    return intens_ind, intens_err_ind, x_center_ind, y_center_ind, \
        a_ind, ellipticity_ind, theta_ind

def last_list(my_list):
    if (type(my_list) is list) or (type(my_list) is numpy.ndarray):
        last = my_list[-1]
    else:
        last = my_list
    return last

def float_dict(dictionary):
    keys = dictionary.keys()
    for n in range(len(keys)):
        for i in range(len(dictionary[keys[n]])):
            try:
                dictionary[keys[n]][i] = float(dictionary[keys[n]][i])
            except ValueError:
                pass
        dictionary[keys[n]] = numpy.asarray(dictionary[keys[n]])
    return dictionary
 
#Changes INDEF to nan       
def indef_fix(var):
    if var == "INDEF":
        var = float('nan')
    else:
        var = float(var)
    return var

#Retrieves a number of variables from an ellipse output table
def get_variables(lines, start_ind):
    num_ellipses = len(lines) - start_ind - 1  #The number of ellipses
    intens, intens_err, x_center, y_center, a, b, theta = \
        init_ellipse_params(num_ellipses)
    indices = ellipse_param_inds(lines)    #finding the indices for the lists
    for n in range(num_ellipses):
        split = lines[n+start_ind].split()
        intens[n] = indef_fix(split[indices[0]])
        intens_err[n] = indef_fix(split[indices[1]])
        x_center[n] = indef_fix(split[indices[2]]) - 1
        y_center[n] = indef_fix(split[indices[3]]) - 1 #iraf counts 1...
        a[n] = indef_fix(split[indices[4]])
        ellipticity = indef_fix(split[indices[5]])
        b[n] = a[n]*(1 - ellipticity)
        theta[n] = indef_fix(split[indices[6]])
    return intens, intens_err, x_center, y_center, a, b, theta

#Reads in the output from IRAF's ellipse
def read_ellipse(filename): 
    lines = loadfile(filename)
    start_ind = header_row(lines)
    intens, intens_err, x_center, y_center, a, b, theta = \
        get_variables(lines, start_ind)
    return intens, intens_err, x_center, y_center, a, b, theta