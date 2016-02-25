import ellipse_funcs
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy
import pickle
from ellipse_params import Galaxy
import pdb

class PlotSurf(object):
    
    def __init__(self, gal):
        self.fig_file = '../figures/'+gal.Galaxy+'.ch1.eps'
        self.set_latex()
        self.set_axes()
        self.plot_surf_bright(gal)
        self.plot_local_scale_length(gal)
        self.plot_asym(gal)
        
    def set_latex(self):
        plt.rc('text', usetex=True)
        plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amssymb}']
        
    def set_axes(self):
        surf_ax, surf_ax_x, surf_ax_y, scale_length_ax, scale_length_ax_x, asym_ax = \
            self.get_axes()
        self.surf_ax = surf_ax
        self.surf_ax_x = surf_ax_x
        self.surf_ax_y = surf_ax_y
        self.scale_length_ax = scale_length_ax
        self.scale_length_ax_x = scale_length_ax_x
        self.asym_ax = asym_ax
    
    def get_axes(self):
        gs1 = gridspec.GridSpec(7, 7) #Determines the spacing
        gs1.update(hspace = 0.3) #Removes space between plots
        surf_ax = plt.subplot(gs1[0:3, :])
        scale_length_ax = plt.subplot(gs1[3:5, :])
        asym_ax = plt.subplot(gs1[5:7, :])
        surf_ax_x = surf_ax.twinx()
        surf_ax_y = surf_ax.twiny()
        scale_length_ax_x = scale_length_ax.twinx()
        return surf_ax, surf_ax_x, surf_ax_y, scale_length_ax, scale_length_ax_x, asym_ax
        
    def plot_surf_bright(self, gal):
        self.plot_surf_bright_main(gal)
        self.plot_surf_bright_x(gal)
        self.plot_surf_bright_y(gal)
        
    def plot_surf_bright_main(self, gal):
        #The main plot
        self.surf_ax.errorbar(gal.arcmin_sma, gal.mag, \
            yerr = [gal.mag_lower, gal.mag_upper], \
            fmt = 'ko', capsize = 0, markersize = 4.0)
        
        #Various lines    
        self.surf_ax.text(1.025*gal.r25, min(gal.mag[1:]) + 0.4, 'R25') 
        self.surf_ax.plot([gal.r25, gal.r25], [0, 100], 'k')
        self.surf_ax.plot([gal.bounds[0], gal.bounds[0]], 
            [0, 100], 'k', ls = 'dashed')
        self.surf_ax.plot([gal.bounds[-1], gal.bounds[-1]], 
            [0, 100], 'k', ls = 'dashed')
        if gal.r_bar != 0:
            self.surf_ax.plot([gal.r_bar, gal.r_bar], [0,100], 'k')
            self.surf_ax.plot([gal.r_bar*2, gal.r_bar*2.0], [0,100], 'k')
        for r_break in gal.r_break:
            self.surf_ax.plot([r_break, r_break], [0, 100])
        
        #Setting text and limits    
        auto_xlim(self.surf_ax, max(gal.arcmin_sma))
        plt.setp([self.surf_ax.get_xticklabels()], visible=False)
        self.surf_ax.set_ylim((min(gal.mag[1:]) - 1, max(gal.mag) + 2))
        self.surf_ax.invert_yaxis()
        type_string = type_string_calc(gal.break_type[0])
        self.surf_ax.text(0.1*max(gal.arcmin_sma), 
            max(gal.mag[1:]) - 1.8, gal.Galaxy)
        self.surf_ax.text(0.1*max(gal.arcmin_sma), 
            max(gal.mag[1:]) - 0.4, 'Type-'+type_string)
        self.surf_ax.set_ylabel(r"${\rm \mu~(Vega~mag/\square '')}$")
        self.surf_ax.grid(True)
                
    def plot_surf_bright_x(self, gal):
        self.surf_ax_x.set_ylim([ellipse_funcs.mag2mass(max(gal.mag) + 2), 
            ellipse_funcs.mag2mass(min(gal.mag[1:]) - 1)])
        self.surf_ax_x.set_ylabel(r'${\rm log_{10}~\Sigma_{\star}~(M_{\odot}~pc^{-2})}$')
        
    def plot_surf_bright_y(self, gal):
        kpc_max = ellipse_funcs.arcmin2kpc(max(gal.arcmin_sma), gal.d_kpc)
        auto_xlim(self.surf_ax_y, kpc_max)
        self.surf_ax_y.set_xlabel(r'${\rm Semi-Major~Axis~(kpc)}$')
        
    def plot_local_scale_length(self, gal):
        self.scale_length_ax.errorbar(gal.h_local_sma, gal.h_local, 
            yerr = gal.h_local_unc, 
            fmt = 'ko', capsize = 0, markersize = 4.0)
        self.scale_length_ax.plot([gal.r25, gal.r25], [0, 100], 'k')
        self.scale_length_ax.plot([0, 100], [gal.scale_med, gal.scale_med], 'k')
        auto_xlim(self.scale_length_ax, max(gal.arcmin_sma))
        ind = numpy.where(gal.h_local_unc/gal.h_local < 1.0)
        y_lim = scale_ylim(gal.h_local[ind], gal.h_local_unc[ind])
        if y_lim > 5.0:
            y_lim = 5.0
        self.scale_length_ax.text(0.85*max(gal.arcmin_sma[ind]), 
            gal.scale_med + 0.05*y_lim, 'Median')
        self.scale_length_ax.set_ylim((0, y_lim))
        self.scale_length_ax.set_ylabel(r"$h_{{\rm local}}~(')$")
        self.scale_length_ax.grid(True)
        self.plot_local_scale_length_x(gal, y_lim)
        remove_xlabels([self.scale_length_ax, self.scale_length_ax_x])

    def plot_local_scale_length_x(self, gal, y_lim):
        kpc_y_lim = ellipse_funcs.arcmin2kpc(y_lim, gal.d_kpc)
        self.scale_length_ax_x.set_ylim((0, kpc_y_lim))
        self.scale_length_ax_x.set_ylabel(r"$h_{{\rm local}}~({\rm kpc})$")
        
    def plot_asym(self, gal):
        self.asym_ax.plot(gal.arcmin_sma, gal.asym, 'ko', markersize = 4.0)
        auto_xlim(self.asym_ax, max(gal.arcmin_sma))
        self.asym_ax.set_xlabel(r"${\rm Semi-Major~Axis~(}${\tt '}${\rm )}$")
        self.asym_ax.set_ylabel(r'$a$')
        self.asym_ax.grid(True)
        
def remove_xlabels(axes):
    for ax in axes:
        plt.setp([ax.get_xticklabels()], visible=False)
        
def auto_xlim(ax, max_x):
    ax.set_xlim((0, 1.1*max_x))
    
def scale_ylim(scale_length, scale_length_unc):
    ind = numpy.where(numpy.isfinite(scale_length + scale_length_unc))
    max_length = 1.1*max(scale_length[ind] + scale_length_unc[ind])
    if max_length - int(max_length) < 0.5:
        y_lim = max_length + 0.5
    else:
        y_lim = max_length
    if max_length > 10:
        y_lim = 10.5
    return y_lim

def type_string_calc(break_type):
    if break_type == 1:
        type_string = 'I'
    if int(break_type) == 2:
        type_string = 'II'
    if break_type == 3:
        type_string = 'III'
    return type_string   
        
def main():
    with open('../pickles/gals.pickle', 'r') as handle:
      gals = pickle.load(handle)
    for gal in gals:
        print gal.Galaxy
        plot = PlotSurf(gal)
        plt.savefig(plot.fig_file)
    
if __name__ == '__main__':
    main()