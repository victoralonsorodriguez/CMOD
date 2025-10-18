
import numpy as np
import sys
import pdb

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import functools

from py_open_fits import open_fits
from py_convert_functions import values_counts_to_mag,px_to_kpc,kpc_to_px


def plot_images(gal_img,
                fig_name,
                cons,
                path=False,
                ellip=False,isolist=None,
                mag=False,res=False,
                conver_to_mag=False,
                counts=True,
                output_path='.',
                log_scale = False):
    
    # Enable LaTeX rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')  # Use a serif font for LaTeX rendering
    plt.rc('font', size=14)  # Adjust size to your preference
    # Define the LaTeX preamble with siunitx
    plt.rcParams['text.latex.preamble'] = r'''
                \usepackage{siunitx}
                \sisetup{
                  detect-family,
                  separate-uncertainty=true,
                  output-decimal-marker={.},
                  exponent-product=\cdot,
                  inter-unit-product=\cdot,
                }
                \DeclareSIUnit{\cts}{cts}
                '''
    
    if path == True:
        
        print('Path was passed to function to plot the image')
        _,gal_img_mag,_ = open_fits(gal_img)
        
    else:
            
        gal_img_mag = gal_img
        
    if counts==True and conver_to_mag==True:
        gal_img_mag = values_counts_to_mag(gal_img_mag,cons)
        
    if counts==True and log_scale == True:
        gal_img_mag = np.log10(gal_img_mag)
    
    x_len = gal_img_mag.shape[1]
    y_len = gal_img_mag.shape[0]
    
    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(1,1,1)
    
    if conver_to_mag==True or mag==True:
        color_map = 'viridis_r'
        bar_label = r'$\mu\ [\mathrm{mag/arcsec^2}]$'
        
    elif counts==True and conver_to_mag==False:
        color_map = 'viridis'
        bar_label = r'$Intensisty\, [\mathrm{counts}]$'
        
        if log_scale == True:
            bar_label = r'$\log_{10}(Intensisty/\mathrm{counts})$'
    
    if res == True:
        im = ax.imshow(gal_img_mag, origin='lower',cmap=color_map,
                       vmin=-1,vmax=1)
    else: 
        im = ax.imshow(gal_img_mag, origin='lower',cmap=color_map)
    plt.colorbar(im, pad=0.11)
    
    if conver_to_mag==True:
        fig.axes[1].invert_yaxis()
    
    fig.axes[1].set(ylabel=bar_label)
    fig.axes[1].minorticks_on()
    
    if ellip == True:

        sma_max = int(max(isolist.to_table()['sma']))
        smas = np.linspace(5, sma_max, 30)
        
        for sma in smas:
            iso = isolist.get_closest(sma)
            x, y, = iso.sampled_coordinates()
            plt.plot(x, y, color='white',linewidth=0.5)
            
    
    # Customizing the plots       
    # X bottom axis is common
    ax.set_xlabel(r'$X\ dimension\ [\mathrm{pix}]$')
    ax.set_xlim(left=0.0)
    
    ax.set_xticks(np.linspace(0, x_len, 7))
    ax.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    
    #ax.grid(True,alpha=0.2)
    
    # TOP X axis is common
    axxtop = ax.secondary_xaxis('top',
                                functions=(functools.partial(px_to_kpc,inst=cons[0]),
                                           functools.partial(kpc_to_px,inst=cons[0])))
    
    px_ticks = ax.get_xticks()
    arcsec_ticks = px_to_kpc(px_ticks,inst=cons[0])
    axxtop.set_xticks(arcsec_ticks)
    
    axxtop.minorticks_on()
    
    axxtop.set_xlabel(r'$X\ dimension\ [\mathrm{kpc}]$',labelpad=8)
    axxtop.tick_params(axis='x', which='major')
    
    # Y left axis
    ax.set_ylabel(r'$Y\ dimension\ [\mathrm{pix}]$')
    ax.set_yticks(np.linspace(0, y_len, 7))
    
    ax.minorticks_on()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_ylim(bottom=0.0)
    
    # Y right axis
    
    axyrig = ax.secondary_yaxis('right',
                                functions=(functools.partial(px_to_kpc,inst=cons[0]),
                                           functools.partial(kpc_to_px,inst=cons[0])))
    
    px_ticks = ax.get_yticks()
    arcsec_ticks = px_to_kpc(px_ticks,inst=cons[0])
    axyrig.set_yticks(arcsec_ticks)
    
    axyrig.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axyrig.yaxis.set_minor_formatter(plt.NullFormatter())
    
    axyrig.minorticks_on()
    
    axyrig.set_ylabel(r'$Y\ dimension\ [\mathrm{kpc}]$')
    axyrig.tick_params(axis='y', which='major')
    
    if '-png' in sys.argv:
        img_format = 'png'
    else: 
        img_format = 'pdf'

    fig_name = f'{fig_name}.{img_format}'
    fig_path = f'{output_path}/{fig_name}'
    plt.savefig(f'{fig_path}',format=img_format, dpi=1000, bbox_inches='tight')
    plt.close()