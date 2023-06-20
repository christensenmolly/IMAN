
#!/usr/bin/python
import numpy as np
import math
import scipy.special   # type: ignore
import scipy.optimize   # type: ignore
from uncertainties import ufloat
from uncertainties.umath import *
import sys
import collections
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from scipy import stats
import warnings
import matplotlib
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
from matplotlib.ticker import MaxNLocator
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')
mpl.rcParams['errorbar.capsize'] = 3
from astropy.stats import sigma_clipped_stats
from matplotlib.pyplot import cm
 
warnings.filterwarnings("ignore")
 
~/MyGit/IMAN/
 
 
sys.path.append('~/MyGit/IMAN/misc_funcs')
sys.path.append('~/MyGit/IMAN/plotting/graphics')
 
import read_data
import graph_plotter
 
from collections import OrderedDict
 
import func_prepare_results
 
 
 
def plot_h_lambda(tex_lines, names, Res_disk_main, h_100, h_160, h_250, h_350, h_500, ouput_file='h_lambda.png'):
        h_100 = np.array(h_100)
        h_160 = np.array(h_160)
        h_250 = np.array(h_250)
        h_350 = np.array(h_350)
        h_500 = np.array(h_500)
    
        h_100_aver,h_100_med,h_100_std = sigma_clipped_stats(h_100/h_100, sigma=2.0)
        h_160_aver,h_160_med,h_160_std = sigma_clipped_stats(h_160/h_100, sigma=2.0)
        h_250_aver,h_250_med,h_250_std = sigma_clipped_stats(h_250/h_100, sigma=2.0)
        h_350_aver,h_350_med,h_350_std = sigma_clipped_stats(h_350/h_100, sigma=2.0)
        h_500_aver,h_500_med,h_500_std = sigma_clipped_stats(h_500/h_100, sigma=2.0)        
 
        print('\t h_100/h100 = $%.2f(%.2f) \pm %.2f$' % (h_100_aver, h_100_med, h_100_std))
        print('\t h_160/h100 = $%.2f(%.2f) \pm %.2f$' % (h_160_aver, h_160_med, h_160_std))
        print('\t h_250/h100 = $%.2f(%.2f) \pm %.2f$' % (h_250_aver, h_250_med, h_250_std))
        print('\t h_350/h100 = $%.2f(%.2f) \pm %.2f$' % (h_350_aver, h_350_med, h_350_std))
        print('\t h_500/h100 = $%.2f(%.2f) \pm %.2f$' % (h_500_aver, h_500_med, h_500_std))
 
 
        tex_lines['100'] = tex_lines['100'] + '& $%.2f(%.2f) \pm %.2f$' % (h_100_aver, h_100_med, h_100_std)
        tex_lines['160'] = tex_lines['160'] + '& $%.2f(%.2f) \pm %.2f$' % (h_160_aver, h_160_med, h_160_std)
        tex_lines['250'] = tex_lines['250'] + '& $%.2f(%.2f) \pm %.2f$' % (h_250_aver, h_250_med, h_250_std)
        tex_lines['350'] = tex_lines['350'] + '& $%.2f(%.2f) \pm %.2f$' % (h_350_aver, h_350_med, h_350_std)
        tex_lines['500'] = tex_lines['500'] + '& $%.2f(%.2f) \pm %.2f$' % (h_500_aver, h_500_med, h_500_std)
 
        colors=iter(cm.nipy_spectral(np.linspace(0,1,len(h_100))))
        outer_colors = collections.OrderedDict()
        
        for k in range(len(h_100)):
            Color = next(colors)
            res = [h_100[k]/h_100[k],h_160[k]/h_100[k],h_250[k]/h_100[k],h_350[k]/h_100[k],h_500[k]/h_100[k]]
 
            wavelengths = np.array([100.,160.,250.,350.,500.], float)
            ax,overlayed_number = graph_plotter.main(wavelengths,res,x_err=None,y_err=None,x_dim=[80.,600.],y_dim=[0.81,2.5],
                    figure=[None,0,None,0,[5,5],False],
                    text=None,
                    regression_line=None,
                    marker=['-',5,Color,Color],
                    line=['-',1,Color],
                    labels=['$\lambda$ [$\mu$m]','$h_{R,\lambda}/h_{R,100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
            if names[k] in ['NGC4217','NGC4302','NGC4437','NGC5907']:
                outer_colors[names[k]] = Color
 
        for k in range(len(Res_disk_main["Name"])):
            if Res_disk_main["Name"][k] in ['NGC4217','NGC4302','NGC4437','NGC5907']:
                res = [Res_disk_main['h2_100'][k]/Res_disk_main['h2_100'][k], Res_disk_main['h2_160'][k]/Res_disk_main['h2_100'][k], Res_disk_main['h2_250'][k]/Res_disk_main['h2_100'][k], Res_disk_main['h2_350'][k]/Res_disk_main['h2_100'][k], Res_disk_main['h2_500'][k]/Res_disk_main['h2_100'][k]]
                ax.plot(wavelengths, res, ls='--', lw=2, color='black')
                ax.plot(wavelengths, res, ls='--', lw=1, color=outer_colors[Res_disk_main["Name"][k]])
 
 
        h_aver = np.array([h_100_aver,h_160_aver,h_250_aver,h_350_aver,h_500_aver])
        h_std  = np.array([h_100_std,h_160_std,h_250_std,h_350_std,h_500_std])
        
        ax.axhline(y=1.0, ls="-.", lw=1, color='black')
 
        ax,overlayed_number = graph_plotter.main(wavelengths,h_aver,x_err=None,y_err=h_std,x_dim=[80.,600.],y_dim=[0.81,2.5],
                    figure=[ax,overlayed_number,ouput_file,0,[5,5],True],
                    text=None,
                    regression_line=None,
                    marker=['-',5,'black','black'],
                    line=['-',3,'black'],
                    labels=['$\lambda$ [$\mu$m]','$h_{R,\lambda}/h_{R,100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
 
def plot_z0_lambda(tex_lines, Res_disk_main, Res_sersic_main, ouput_file='z0_lambda.png'):
        '''
        for k in range(len(names)):
            if not func_prepare_results.check_goodness(round(h_100[k]/Scale[k],2),'100'):
                h_100[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_160[k]/Scale[k],2),'160'):
                h_160[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_250[k]/Scale[k],2),'250'):
                h_250[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_350[k]/Scale[k],2),'350'):
                h_350[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_500[k]/Scale[k],2),'500'):
                h_500[k] = float('nan')
        '''
        
        names = Res_sersic_main["Name"]
 
        h_100 = Res_disk_main["z0_100"]
        h_160 = Res_disk_main["z0_160"]
        h_250 = Res_disk_main["z0_250"]
        h_350 = Res_disk_main["z0_350"]
        h_500 = Res_disk_main["z0_500"]
 
        for k in range(len(names)):
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_100'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_100'][k], Res_sersic_main['q_100'][k],'100'):
                h_100[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_160'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_160'][k], Res_sersic_main['q_160'][k],'160'):
                h_160[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_250'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_250'][k], Res_sersic_main['q_250'][k],'250'):
                h_250[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_350'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_350'][k], Res_sersic_main['q_350'][k],'350'):
                h_350[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_500'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_500'][k], Res_sersic_main['q_500'][k],'500'):
                h_500[k] = float('nan')
 
 
                
                
        names = np.array(names)
        h_100 = np.array(h_100)
        h_160 = np.array(h_160)
        h_250 = np.array(h_250)
        h_350 = np.array(h_350)
        h_500 = np.array(h_500)
        
        
        
        h_100_aver,h_100_med,h_100_std = sigma_clipped_stats(h_100/h_100, sigma=3.0)
        h_160_aver,h_160_med,h_160_std = sigma_clipped_stats(h_160/h_100, sigma=3.0)
        h_250_aver,h_250_med,h_250_std = sigma_clipped_stats(h_250/h_100, sigma=3.0)
        h_350_aver,h_350_med,h_350_std = sigma_clipped_stats(h_350/h_100, sigma=3.0)
        h_500_aver,h_500_med,h_500_std = sigma_clipped_stats(h_500/h_100, sigma=3.0)        
 
        print('\t z0_100/z0_100 = $%.2f(%.2f) \pm %.2f$' % (h_100_aver, h_100_med, h_100_std))
        print('\t z0_160/z0_100 = $%.2f(%.2f) \pm %.2f$' % (h_160_aver, h_160_med, h_160_std))
        print('\t z0_250/z0_h100 = $%.2f(%.2f) \pm %.2f$' % (h_250_aver, h_250_med, h_250_std))
        print('\t z0_350/z0_h100 = $%.2f(%.2f) \pm %.2f$' % (h_350_aver, h_350_med, h_350_std))
        print('\t z0_500/z0_h100 = $%.2f(%.2f) \pm %.2f$' % (h_500_aver, h_500_med, h_500_std))
 
        tex_lines['100'] = tex_lines['100'] + '& $%.2f(%.2f) \pm %.2f$ \\\\ \n' % (h_100_aver, h_100_med, h_100_std)
        tex_lines['160'] = tex_lines['160'] + '& $%.2f(%.2f) \pm %.2f$ \\\\ \n' % (h_160_aver, h_160_med, h_160_std)
        tex_lines['250'] = tex_lines['250'] + '& $%.2f(%.2f) \pm %.2f$ \\\\ \n' % (h_250_aver, h_250_med, h_250_std)
        tex_lines['350'] = tex_lines['350'] + '& $%.2f(%.2f) \pm %.2f$ \\\\ \n' % (h_350_aver, h_350_med, h_350_std)
        tex_lines['500'] = tex_lines['500'] + '& $%.2f(%.2f) \pm %.2f$ \\\\ \n' % (h_500_aver, h_500_med, h_500_std)
    
 
 
        colors=iter(cm.nipy_spectral(np.linspace(0,1,len(h_100))))
 
        for k in range(len(h_100)):
            Color = next(colors)
            res = [h_100[k]/h_100[k],h_160[k]/h_100[k],h_250[k]/h_100[k],h_350[k]/h_100[k],h_500[k]/h_100[k]]
 
            wavelengths = np.array([100.,160.,250.,350.,500.], float)
            ax,overlayed_number = graph_plotter.main(wavelengths,res,x_err=None,y_err=None,x_dim=[80.,600.],y_dim=[0.81,2.49],
                    figure=[None,0,None,0,[5,5],False],
                    text=None,
                    regression_line=None,
                    marker=['-',5,Color,Color],
                    line=['-',1,Color],
                    labels=['$\lambda$ [$\mu$m]','$h_{z,\lambda}/h_{z,100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
 
 
        h_aver = np.array([h_100_aver,h_160_aver,h_250_aver,h_350_aver,h_500_aver])
        h_std  = np.array([h_100_std,h_160_std,h_250_std,h_350_std,h_500_std])
        
        ax.axhline(y=1.0, ls="-.", lw=1, color='black')
 
        ax,overlayed_number = graph_plotter.main(wavelengths,h_aver,x_err=None,y_err=h_std,x_dim=[80.,600.],y_dim=[0.81,2.49],
                    figure=[ax,overlayed_number,ouput_file,0,[5,5],True],
                    text=None,
                    regression_line=None,
                    marker=['-',5,'black','black'],
                    line=['-',3,'black'],
                    labels=['$\lambda$ [$\mu$m]','$h_{z,\lambda}/h_{z,100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
        
        
def plot_re_lambda(tex_lines, names, re_100, re_160, re_250, re_350, re_500, ouput_file='re_lambda.png'):
        re_100 = np.array(re_100)
        re_160 = np.array(re_160)
        re_250 = np.array(re_250)
        re_350 = np.array(re_350)
        re_500 = np.array(re_500)
    
        re_100_aver, re_100_med, re_100_std = sigma_clipped_stats(re_100 / re_100, sigma=3.0)
        re_160_aver, re_160_med, re_160_std = sigma_clipped_stats(re_160 / re_100, sigma=3.0)
        re_250_aver, re_250_med, re_250_std = sigma_clipped_stats(re_250 / re_100, sigma=3.0)
        re_350_aver, re_350_med, re_350_std = sigma_clipped_stats(re_350 / re_100, sigma=3.0)
        re_500_aver, re_500_med, re_500_std = sigma_clipped_stats(re_500 / re_100, sigma=3.0)
 
        print('IMFIT (Usachev):' )
        print('\t re_100 = $%.2f(%.2f) \pm %.2f$' % (re_100_aver, re_100_med, re_100_std))
        print('\t re_160 = $%.2f(%.2f) \pm %.2f$' % (re_160_aver, re_160_med, re_160_std)) 
        print('\t re_250 = $%.2f(%.2f) \pm %.2f$' % (re_250_aver, re_250_med, re_250_std))
        print('\t re_350 = $%.2f(%.2f) \pm %.2f$' % (re_350_aver, re_350_med, re_350_std))
        print('\t re_500 = $%.2f(%.2f) \pm %.2f$\n\n' % (re_500_aver, re_500_med, re_500_std))    
 
        tex_lines['100'] = tex_lines['100'] + '& $%.2f(%.2f) \pm %.2f$' % (re_100_aver, re_100_med, re_100_std)
        tex_lines['160'] = tex_lines['160'] + '& $%.2f(%.2f) \pm %.2f$' % (re_160_aver, re_160_med, re_160_std)
        tex_lines['250'] = tex_lines['250'] + '& $%.2f(%.2f) \pm %.2f$' % (re_250_aver, re_250_med, re_250_std)
        tex_lines['350'] = tex_lines['350'] + '& $%.2f(%.2f) \pm %.2f$' % (re_350_aver, re_350_med, re_350_std)
        tex_lines['500'] = tex_lines['500'] + '& $%.2f(%.2f) \pm %.2f$' % (re_500_aver, re_500_med, re_500_std)
    
        colors=iter(cm.nipy_spectral(np.linspace(0,1,len(re_100))))
 
        for k in range(len(re_100)):
            Color = next(colors)
            res = [re_100[k]/re_100[k],
                re_160[k]/re_100[k],
                re_250[k]/re_100[k],
                re_350[k]/re_100[k],
                re_500[k]/re_100[k]]
 
            wavelengths = np.array([100.,160.,250.,350.,500.], float)
            ax,overlayed_number = graph_plotter.main(wavelengths,res,x_err=None,y_err=None,x_dim=[80.,600.],y_dim=[0.81,2.8],
                    figure=[None,0,None,0,[5,5],False],
                    text=None,
                    regression_line=None,
                    marker=['-',5,Color,Color],
                    line=['-',1,Color],
                    labels=['$\lambda$ [$\mu$m]','$r_{\mathrm{e},\lambda}/r_{\mathrm{e},100}$',18,'black'],
                    legend=[names[k],2, 8, 'white'],
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
        #[np.array(results_pavel['galaxyName'],str)[P_100[k]],2, 6, 'white']
 
        re_aver = [re_100_aver,re_160_aver,re_250_aver,re_350_aver,re_500_aver]
        re_std = [re_100_std,re_160_std,re_250_std,re_350_std,re_500_std]
        
        ax.axhline(y=1.0, ls="-.", lw=1, color='black')
 
        ax,overlayed_number = graph_plotter.main(wavelengths,re_aver,x_err=None,y_err=re_std,x_dim=[80.,600.],y_dim=[0.81,2.8],
                    figure=[ax,overlayed_number,ouput_file,0,[5,5],True],
                    text=None,
                    regression_line=None,
                    marker=['-',5,'black','black'],
                    line=['-',3,'black'],
                    labels=['$\lambda$ [$\mu$m]','$r_{\mathrm{e},\lambda}/r_{\mathrm{e},100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
 
 
 
 
def plot_q_lambda(tex_lines, Res_disk_main, Res_sersic_main, ouput_file='q_lambda.png'):
        names = Res_sersic_main["Name"]
    
        '''
        for k in range(len(names)):
            if not func_prepare_results.check_goodness(round(hz_100[k]/Scale[k],2),'100'):
                q_100[k] = float('nan')
            if not func_prepare_results.check_goodness(round(hz_160[k]/Scale[k],2),'160'):
                q_160[k] = float('nan')
            if not func_prepare_results.check_goodness(round(hz_250[k]/Scale[k],2),'250'):
                q_250[k] = float('nan')
            if not func_prepare_results.check_goodness(round(hz_350[k]/Scale[k],2),'350'):
                q_350[k] = float('nan')
            if not func_prepare_results.check_goodness(round(hz_500[k]/Scale[k],2),'500'):
                q_500[k] = float('nan')
        '''
        
        q_100 = Res_sersic_main["q_100"]
        q_160 = Res_sersic_main["q_160"]
        q_250 = Res_sersic_main["q_250"]
        q_350 = Res_sersic_main["q_350"]
        q_500 = Res_sersic_main["q_500"]
 
        for k in range(len(names)):
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_100'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_100'][k], Res_sersic_main['q_100'][k],'100'):
                q_100[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_160'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_160'][k], Res_sersic_main['q_160'][k],'160'):
                q_160[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_250'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_250'][k], Res_sersic_main['q_250'][k],'250'):
                q_250[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_350'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_350'][k], Res_sersic_main['q_350'][k],'350'):
                q_350[k] = float('nan')
            if not func_prepare_results.check_goodness_ser(Res_disk_main['Name'][k], round(Res_sersic_main['re_500'][k]/Res_sersic_main['Scale'][k],2), Res_sersic_main['n_500'][k], Res_sersic_main['q_500'][k],'500'):
                q_500[k] = float('nan')
 
 
 
        q_100 = np.array(q_100)
        q_160 = np.array(q_160)
        q_250 = np.array(q_250)
        q_350 = np.array(q_350)
        q_500 = np.array(q_500)
 
        q_100_aver,q_100_med,q_100_std = sigma_clipped_stats(q_100/q_100, sigma=2.0)
 
        q_160_aver,q_160_med,q_160_std = sigma_clipped_stats(q_160/q_100, sigma=2.0)
 
        q_250_aver,q_250_med,q_250_std = sigma_clipped_stats(q_250/q_100, sigma=2.0)
 
        q_350_aver,q_350_med,q_350_std = sigma_clipped_stats(q_350/q_100, sigma=2.0)
 
        q_500_aver,q_500_med,q_500_std = sigma_clipped_stats(q_500/q_100, sigma=2.0)
 
 
        print('\t q_100/q_100= = $%.2f(%.2f) \pm %.2f$' % (q_100_aver, q_100_med, q_100_std))
        print('\t q_160/q_100= = $%.2f(%.2f) \pm %.2f$' % (q_160_aver, q_160_med, q_160_std))
        print('\t q_250/q_100= = $%.2f(%.2f) \pm %.2f$' % (q_250_aver, q_250_med, q_250_std))
        print('\t q_350/q_100= = $%.2f(%.2f) \pm %.2f$' % (q_350_aver, q_350_med, q_350_std))
        print('\t q_500/q_100= = $%.2f(%.2f) \pm %.2f$\n\n' % (q_500_aver, q_500_med, q_500_std))    
        
        tex_lines['100'] = tex_lines['100'] + '& $%.2f(%.2f) \pm %.2f$' % (q_100_aver, q_100_med, q_100_std)
        tex_lines['160'] = tex_lines['160'] + '& $%.2f(%.2f) \pm %.2f$' % (q_160_aver, q_160_med, q_160_std)
        tex_lines['250'] = tex_lines['250'] + '& $%.2f(%.2f) \pm %.2f$' % (q_250_aver, q_250_med, q_250_std)
        tex_lines['350'] = tex_lines['350'] + '& $%.2f(%.2f) \pm %.2f$' % (q_350_aver, q_350_med, q_350_std)
        tex_lines['500'] = tex_lines['500'] + '& $%.2f(%.2f) \pm %.2f$' % (q_500_aver, q_500_med, q_500_std)
 
        colors=iter(cm.nipy_spectral(np.linspace(0,1,len(q_100))))
 
        for k in range(len(q_100)):
            Color = next(colors)
            res = [q_100[k]/q_100[k],
                q_160[k]/q_100[k],
                q_250[k]/q_100[k],
                q_350[k]/q_100[k],
                q_500[k]/q_100[k]]
            #print(names[k], q_500[k],q_100[k])
            wavelengths = np.array([100.,160.,250.,350.,500.], float)
            ax,overlayed_number = graph_plotter.main(wavelengths,res,x_err=None,y_err=None,x_dim=[80.,600.],y_dim=[0.81,1.39],
                    figure=[None,0,None,0,[5,5],False],
                    text=None,
                    regression_line=None,
                    marker=['-',5,Color,Color],
                    line=['-',1,Color],
                    labels=['$\lambda$ [$\mu$m]','$q_{\lambda}/q_{100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
 
 
        q_aver = np.array([q_100_aver,q_160_aver,q_250_aver,q_350_aver,q_500_aver])
        q_std  = np.array([q_100_std,q_160_std,q_250_std,q_350_std,q_500_std])
        
        ax.axhline(y=1.0, ls="-.", lw=1, color='black')
 
        ax,overlayed_number = graph_plotter.main(wavelengths,q_aver,x_err=None,y_err=q_std,x_dim=[80.,600.],y_dim=[0.81,1.39],
                    figure=[ax,overlayed_number,ouput_file,0,[5,5],True],
                    text=None,
                    regression_line=None,
                    marker=['-',5,'black','black'],
                    line=['-',3,'black'],
                    labels=['$\lambda$ [$\mu$m]','$q_{\lambda}/q_{100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
 
 
 
 
 
 
def plot_z0h_lambda(names, Scale, h_100, h_160, h_250, h_350, h_500, z0_100, z0_160, z0_250, z0_350, z0_500, ouput_file='z0h_lambda.png'):
        q_100 = np.array(z0_100)/np.array(h_100)
        q_160 = np.array(z0_160)/np.array(h_160)
        q_250 = np.array(z0_250)/np.array(h_250)        
        q_350 = np.array(z0_350)/np.array(h_350)        
        q_500 = np.array(z0_500)/np.array(h_500)        
        
        
        
        for k in range(len(names)):
            if not func_prepare_results.check_goodness(round(h_100[k]/Scale[k],2),'100'):
                q_100[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_160[k]/Scale[k],2),'160'):
                q_160[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_250[k]/Scale[k],2),'250'):
                q_250[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_350[k]/Scale[k],2),'350'):
                q_350[k] = float('nan')
            if not func_prepare_results.check_goodness(round(h_500[k]/Scale[k],2),'500'):
                q_500[k] = float('nan')
 
 
        q_100 = np.array(q_100)
        q_160 = np.array(q_160)
        q_250 = np.array(q_250)
        q_350 = np.array(q_350)
        q_500 = np.array(q_500)
 
        q_100_aver,q_100_med,q_100_std = sigma_clipped_stats(q_100/q_100, sigma=2.0)
 
        q_160_aver,q_160_med,q_160_std = sigma_clipped_stats(q_160/q_100, sigma=2.0)
 
        q_250_aver,q_250_med,q_250_std = sigma_clipped_stats(q_250/q_100, sigma=2.0)
 
        q_350_aver,q_350_med,q_350_std = sigma_clipped_stats(q_350/q_100, sigma=2.0)
 
        q_500_aver,q_500_med,q_500_std = sigma_clipped_stats(q_500/q_100, sigma=2.0)
 
 
        print('\t h/z0_100/h/z0_100= = $%.2f(%.2f) \pm %.2f$' % (q_100_aver, q_100_med, q_100_std))
        print('\t h/z0_160/h/z0_100= = $%.2f(%.2f) \pm %.2f$' % (q_160_aver, q_160_med, q_160_std))
        print('\t h/z0_250/h/z0_100= = $%.2f(%.2f) \pm %.2f$' % (q_250_aver, q_250_med, q_250_std))
        print('\t h/z0_350/h/z0_100= = $%.2f(%.2f) \pm %.2f$' % (q_350_aver, q_350_med, q_350_std))
        print('\t h/z0_500/h/z0_100= = $%.2f(%.2f) \pm %.2f$\n\n' % (q_500_aver, q_500_med, q_500_std))    
 
 
        colors=iter(cm.nipy_spectral(np.linspace(0,1,len(q_100))))
 
        for k in range(len(q_100)):
            Color = next(colors)
            res = [q_100[k]/q_100[k],
                q_160[k]/q_100[k],
                q_250[k]/q_100[k],
                q_350[k]/q_100[k],
                q_500[k]/q_100[k]]
 
            wavelengths = np.array([100.,160.,250.,350.,500.], float)
            ax,overlayed_number = graph_plotter.main(wavelengths,res,x_err=None,y_err=None,x_dim=[80.,600.],y_dim=[0.51,1.39],
                    figure=[None,0,None,0,[5,5],False],
                    text=None,
                    regression_line=None,
                    marker=['-',5,Color,Color],
                    line=['-',1,Color],
                    labels=['$\lambda$ [$\mu$m]','$h_{z,\lambda}/h_{R,\lambda}\, / \,h_{z,100}/h_{R,100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)
 
 
 
        q_aver = np.array([q_100_aver,q_160_aver,q_250_aver,q_350_aver,q_500_aver])
        q_std  = np.array([q_100_std,q_160_std,q_250_std,q_350_std,q_500_std])
        
        ax.axhline(y=1.0, ls="-.", lw=1, color='black')
 
        ax,overlayed_number = graph_plotter.main(wavelengths,q_aver,x_err=None,y_err=q_std,x_dim=[80.,600.],y_dim=[0.51,1.39],
                    figure=[ax,overlayed_number,ouput_file,0,[5,5],True],
                    text=None,
                    regression_line=None,
                    marker=['-',5,'black','black'],
                    line=['-',3,'black'],
                    labels=['$\lambda$ [$\mu$m]','$h_{z,\lambda}/h_{R,\lambda}\, / \,h_{z,100}/h_{R,100}$',18,'black'],
                    legend=None,
                    diag=False,
                    xticks=np.array(wavelengths, int),
                    yticks=None,
                    add_letter=None,
                    add_line=None,
                    add_points=None,
                    x_logscale=True,
                    y_logscale=False,
                    capsize=4.)



results,units = read_data.main('results_all.txt', delimiter='\t', header_line=0, units_line=None, skip_lines = [], comment=None)

galaxies = results['galaxy']
dust = results['dust']


re_norm_galaxies = []

for k in range(len(galaxies)):
    
        


