#!/usr/bin/python
import sys
import collections
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
import numpy as np
from math import *
from scipy import stats
import warnings
import matplotlib
import matplotlib.gridspec as gridspec
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

warnings.filterwarnings("ignore")




def plot_graphs_I(pars, sel_galaxies=, figsize=(8, 10), color='black', size=4, plot_func=0):
        pars =
        
        
        
        
        plt.figure(0, figsize=figsize)
        
        
        
        gs = gridspec.GridSpec(5, 3)
        gs.update(hspace=0.0, wspace=0.35)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax3 = plt.subplot(gs[2])
        ax4 = plt.subplot(gs[3])
        ax5 = plt.subplot(gs[4])
        ax6 = plt.subplot(gs[5])
        ax7 = plt.subplot(gs[6])
        ax8 = plt.subplot(gs[7])
        ax9 = plt.subplot(gs[8])
        ax10 = plt.subplot(gs[9])
        ax11 = plt.subplot(gs[10])
        ax12 = plt.subplot(gs[11])
        ax13 = plt.subplot(gs[12])
        ax14 = plt.subplot(gs[13])
        ax15 = plt.subplot(gs[14])


        results_w1["r_e"] = np.array(results_w1["r_e"]) * 1.375

        # 100
        pars_w1 = ["r_e","n","q"]
        ax = [ax1,ax2,ax3]
        pars_her = ["Eff_rad_100","Ser_ind_100","q_100"]
        Units = ['arcsec',units_w1["n"],units_w1["q"]]

        x_dims = [[0,100],[0,6],[0,1]]
        y_dims = [[0,100],[0,6],[0,1]]
        for k in range(len(pars_w1)):
        par_w1 = pars_w1[k]
        par_her = pars_her[k]
        AX = ax[k]
        if '_' in par_w1:
            Par = par_w1.split('_')[0] + '$_{'+par_w1.split('_')[-1] + '}$'
        else:
            Par = par_w1   
        label_x = Par + r"$^{W1}$"
        label_y = Par + r"$^{100}$"

        if Units[k]!='':
            label_x = label_x + ' (' + Units[k] + ')'
            label_y = label_y + ' (' + Units[k] + ')'
            
        x_dim = x_dims[k]
        y_dim = y_dims[k]
        graph_plotter.main(results_w1[par_w1],results_her[par_her],x_err=None,y_err=None,x_dim=x_dim,y_dim=y_dim,
                figure=[AX,k,'./graphs/'+par_w1+'_W1_'+par_w1+'_100'+'.png',0,[5,5],False],
                text=None,
                regression_line=['-', 2, 'black'],
                marker=['o',3,'gray','gray'],
                line=['-',2,'blue'],
                labels=[None,label_y,13,'black'],
                legend=None,
                diag=True)


        # 160
        pars_w1 = ["r_e","n","q"]
        ax = [ax4,ax5,ax6]
        pars_her = ["Eff_rad_160","Ser_ind_160","q_160"]
        Units = ['arcsec',units_w1["n"],units_w1["q"]]

        x_dims = [[0,100],[0,6],[0,1]]
        y_dims = [[0,100],[0,6],[0,1]]
        for k in range(len(pars_w1)):
        par_w1 = pars_w1[k]
        par_her = pars_her[k]
        AX = ax[k]
        if '_' in par_w1:
            Par = par_w1.split('_')[0] + '$_{'+par_w1.split('_')[-1] + '}$'
        else:
            Par = par_w1   
        label_x = Par + r"$^{W1}$"
        label_y = Par + r"$^{160}$"

        if Units[k]!='':
            label_x = label_x + ' (' + Units[k] + ')'
            label_y = label_y + ' (' + Units[k] + ')'
            
        x_dim = x_dims[k]
        y_dim = y_dims[k]
        graph_plotter.main(results_w1[par_w1],results_her[par_her],x_err=None,y_err=None,x_dim=x_dim,y_dim=y_dim,
                figure=[AX,k,'./graphs/'+par_w1+'_W1_'+par_w1+'_100'+'.png',0,[5,5],False],
                text=None,
                regression_line=['-', 2, 'black'],
                marker=['o',3,'gray','gray'],
                line=['-',2,'blue'],
                labels=[None,label_y,13,'black'],
                legend=None,
                diag=True)

        # 250
        pars_w1 = ["r_e","n","q"]
        ax = [ax7,ax8,ax9]
        pars_her = ["Eff_rad_250","Ser_ind_250","q_250"]
        Units = ['arcsec',units_w1["n"],units_w1["q"]]

        x_dims = [[0,100],[0,6],[0,1]]
        y_dims = [[0,100],[0,6],[0,1]]
        for k in range(len(pars_w1)):
        par_w1 = pars_w1[k]
        par_her = pars_her[k]
        AX = ax[k]
        if '_' in par_w1:
            Par = par_w1.split('_')[0] + '$_{'+par_w1.split('_')[-1] + '}$'
        else:
            Par = par_w1   
        label_x = Par + r"$^{W1}$"
        label_y = Par + r"$^{250}$"

        if Units[k]!='':
            label_x = label_x + ' (' + Units[k] + ')'
            label_y = label_y + ' (' + Units[k] + ')'
            
        x_dim = x_dims[k]
        y_dim = y_dims[k]
        graph_plotter.main(results_w1[par_w1],results_her[par_her],x_err=None,y_err=None,x_dim=x_dim,y_dim=y_dim,
                figure=[AX,k,'./graphs/'+par_w1+'_W1_'+par_w1+'_100'+'.png',0,[5,5],False],
                text=None,
                regression_line=['-', 2, 'black'],
                marker=['o',3,'gray','gray'],
                line=['-',2,'blue'],
                labels=[None,label_y,13,'black'],
                legend=None,
                diag=True)


        # 350
        pars_w1 = ["r_e","n","q"]
        ax = [ax10,ax11,ax12]
        pars_her = ["Eff_rad_350","Ser_ind_350","q_350"]
        Units = ['arcsec',units_w1["n"],units_w1["q"]]

        x_dims = [[0,100],[0,6],[0,1]]
        y_dims = [[0,100],[0,6],[0,1]]
        for k in range(len(pars_w1)):
        par_w1 = pars_w1[k]
        par_her = pars_her[k]
        AX = ax[k]
        if '_' in par_w1:
            Par = par_w1.split('_')[0] + '$_{'+par_w1.split('_')[-1] + '}$'
        else:
            Par = par_w1   
        label_x = Par + r"$^{W1}$"
        label_y = Par + r"$^{350}$"

        if Units[k]!='':
            label_x = label_x + ' (' + Units[k] + ')'
            label_y = label_y + ' (' + Units[k] + ')'
            
        x_dim = x_dims[k]
        y_dim = y_dims[k]
        graph_plotter.main(results_w1[par_w1],results_her[par_her],x_err=None,y_err=None,x_dim=x_dim,y_dim=y_dim,
                figure=[AX,k,'./graphs/'+par_w1+'_W1_'+par_w1+'_100'+'.png',0,[5,5],False],
                text=None,
                regression_line=['-', 2, 'black'],
                marker=['o',3,'gray','gray'],
                line=['-',2,'blue'],
                labels=[None,label_y,13,'black'],
                legend=None,
                diag=True)



        # 500
        pars_w1 = ["r_e","n","q"]
        ax = [ax13,ax14,ax15]
        pars_her = ["Eff_rad_500","Ser_ind_500","q_500"]
        Units = ['arcsec',units_w1["n"],units_w1["q"]]

        x_dims = [[0,100],[0,6],[0,1]]
        y_dims = [[0,100],[0,6],[0,1]]
        for k in range(len(pars_w1)):
        par_w1 = pars_w1[k]
        par_her = pars_her[k]
        AX = ax[k]
        if '_' in par_w1:
            Par = par_w1.split('_')[0] + '$_{'+par_w1.split('_')[-1] + '}$'
        else:
            Par = par_w1   
        label_x = Par + r"$^{W1}$"
        label_y = Par + r"$^{500}$"

        if Units[k]!='':
            label_x = label_x + ' (' + Units[k] + ')'
            label_y = label_y + ' (' + Units[k] + ')'
            
        x_dim = x_dims[k]
        y_dim = y_dims[k]
        graph_plotter.main(results_w1[par_w1],results_her[par_her],x_err=None,y_err=None,x_dim=x_dim,y_dim=y_dim,
                figure=[AX,k,'./graphs/'+par_w1+'_W1_'+par_w1+'_100'+'.png',0,[5,5],False],
                text=None,
                regression_line=['-', 2, 'black'],
                marker=['o',3,'gray','gray'],
                line=['-',2,'blue'],
                labels=[label_x,label_y,13,'black'],
                legend=None,
                diag=True)


        good_galaxies = []
        for k in range(len(model)):
        if model[k]=='bulge+exp_disc':# or model[k]=='bulge':
            good_galaxies.append(k)

        ax1.plot(np.array(results_w1["r_e"])[good_galaxies], np.array(results_her["Eff_rad_100"])[good_galaxies], 'o', color='red', markersize=3)
        ax2.plot(np.array(results_w1["n"])[good_galaxies], np.array(results_her["Ser_ind_100"])[good_galaxies], 'o', color='red', markersize=3)
        ax3.plot(np.array(results_w1["q"])[good_galaxies], np.array(results_her["q_100"])[good_galaxies], 'o', color='red', markersize=3, zorder = 3)

        ax4.plot(np.array(results_w1["r_e"])[good_galaxies], np.array(results_her["Eff_rad_160"])[good_galaxies], 'o', color='red', markersize=3)
        ax5.plot(np.array(results_w1["n"])[good_galaxies], np.array(results_her["Ser_ind_160"])[good_galaxies], 'o', color='red', markersize=3)
        ax6.plot(np.array(results_w1["q"])[good_galaxies], np.array(results_her["q_160"])[good_galaxies], 'o', color='red', markersize=3, zorder = 3)

        ax7.plot(np.array(results_w1["r_e"])[good_galaxies], np.array(results_her["Eff_rad_250"])[good_galaxies], 'o', color='red', markersize=3)
        ax8.plot(np.array(results_w1["n"])[good_galaxies], np.array(results_her["Ser_ind_250"])[good_galaxies], 'o', color='red', markersize=3)
        ax9.plot(np.array(results_w1["q"])[good_galaxies], np.array(results_her["q_250"])[good_galaxies], 'o', color='red', markersize=3, zorder = 3)

        ax10.plot(np.array(results_w1["r_e"])[good_galaxies], np.array(results_her["Eff_rad_350"])[good_galaxies], 'o', color='red', markersize=3)
        ax11.plot(np.array(results_w1["n"])[good_galaxies], np.array(results_her["Ser_ind_350"])[good_galaxies], 'o', color='red', markersize=3)
        ax12.plot(np.array(results_w1["q"])[good_galaxies], np.array(results_her["q_350"])[good_galaxies], 'o', color='red', markersize=3, zorder = 3)

        ax13.plot(np.array(results_w1["r_e"])[good_galaxies], np.array(results_her["Eff_rad_500"])[good_galaxies], 'o', color='red', markersize=3)
        ax14.plot(np.array(results_w1["n"])[good_galaxies], np.array(results_her["Ser_ind_500"])[good_galaxies], 'o', color='red', markersize=3)
        ax15.plot(np.array(results_w1["q"])[good_galaxies], np.array(results_her["q_500"])[good_galaxies], 'o', color='red', markersize=3, zorder = 3)


        ax1.set_xlim(1,100)
        ax2.set_xlim(0.1,6)
        ax3.set_xlim(0.01,1)

        ax4.set_xlim(1,100)
        ax5.set_xlim(0.1,6)
        ax6.set_xlim(0.01,1)

        ax7.set_xlim(1,100)
        ax8.set_xlim(0.1,6)
        ax9.set_xlim(0.01,1)

        ax10.set_xlim(1,100)
        ax11.set_xlim(0.1,6)
        ax12.set_xlim(0.01,1)

        ax13.set_xlim(1,100)
        ax14.set_xlim(0.1,6)
        ax15.set_xlim(0.01,1)


        ax1.set_ylim(0,99)
        ax2.set_ylim(0,5.9)
        ax3.set_ylim(0,0.99)

        ax4.set_ylim(0,99)
        ax5.set_ylim(0,5.9)
        ax6.set_ylim(0,0.99)

        ax7.set_ylim(0,99)
        ax8.set_ylim(0,5.9)
        ax9.set_ylim(0,0.99)

        ax10.set_ylim(0,99)
        ax11.set_ylim(0,5.9)
        ax12.set_ylim(0,0.99)

        ax13.set_ylim(0,99)
        ax14.set_ylim(0,5.9)
        ax15.set_ylim(0,0.99)


        ax1.tick_params(labelbottom='off')
        ax2.tick_params(labelbottom='off')   
        ax3.tick_params(labelbottom='off')   
        ax4.tick_params(labelbottom='off')
        ax5.tick_params(labelbottom='off')
        ax6.tick_params(labelbottom='off')   
        ax7.tick_params(labelbottom='off')   
        ax8.tick_params(labelbottom='off')
        ax9.tick_params(labelbottom='off')
        ax10.tick_params(labelbottom='off')   
        ax11.tick_params(labelbottom='off')   
        ax12.tick_params(labelbottom='off')


        ax1.tick_params(direction='in')
        ax2.tick_params(direction='in')
        ax3.tick_params(direction='in')
        ax4.tick_params(direction='in')
        ax5.tick_params(direction='in')
        ax6.tick_params(direction='in')
        ax7.tick_params(direction='in')
        ax8.tick_params(direction='in')
        ax9.tick_params(direction='in')
        ax10.tick_params(direction='in')
        ax11.tick_params(direction='in')
        ax12.tick_params(direction='in')


        plt.savefig('./graphs/pars.png', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)   
        plt.clf()
        plt.close() 
