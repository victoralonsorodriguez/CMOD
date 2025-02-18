
import sys
import numpy as np
import functools

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from py_px_kpc import px_to_kpc,kpc_to_px,arcsec_to_kpc,kpc_to_arcsec
from py_rad_deg import deg_to_rad,rad_to_deg
from py_mag_counts_convert import fits_mag_to_counts, fits_counts_to_mag, values_counts_to_mag,values_mag_to_counts
from py_convert_ell_axrat import axrat_to_ell,ell_to_axrat

def plot_initial_conditions(fig_name,
                            x_data_tot,
                            y_data_tot,
                            y_data_func = [],
                            labels = [],
                            const=None,
                            galaxy_folder_path='.'):
    
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
    
    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(1,1,1)
    
    colors = ['gold','deepskyblue','lime','deeppink']
    line_styles = ['-','--','-.','--.']
    
    plt.scatter(x_data_tot,y_data_tot,label='Intensity profile',
                marker='o',s=10,
                linewidth=0.15,edgecolor='black',
                color='black',zorder=1)
    
    for data in range(len(y_data_func)):
        plt.plot(x_data_tot,y_data_func[data],label=labels[data],
            linewidth=1,color=colors[data],
            linestyle=line_styles[data],zorder=3)
        
    plt.legend(loc='upper right')
    
    # Customizing the plots       
    # X bottom axis is common
    ax.set_xlabel(r'$X\ dimension\ [\mathrm{pix}]$')
    ax.set_xlim(left=-5)
    
    ax.set_xticks(np.linspace(0, np.max(x_data_tot), 7))
    ax.minorticks_on()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    
    ax.grid(True,alpha=0.2)
    
    # TOP X axis is common
    axxtop = ax.secondary_xaxis('top',
                                functions=(functools.partial(px_to_kpc,inst=const[0]),
                                           functools.partial(kpc_to_px,inst=const[0])))
    
    px_ticks = ax.get_xticks()
    arcsec_ticks = px_to_kpc(px_ticks,inst=const[0])
    axxtop.set_xticks(arcsec_ticks)
    
    axxtop.minorticks_on()
    
    axxtop.set_xlabel(r'$X\ dimension\ [\mathrm{kpc}]$',labelpad=8)
    axxtop.tick_params(axis='x', which='major')
    
    # Y left axis
    ax.set_ylabel(r'$\ln(Intensity/\mathrm{counts})$')
    ax.set_yticks(np.linspace(np.min(y_data_tot), np.max(y_data_tot), 7))
    
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.minorticks_on()
    
    # Y right axis
    axyrig = ax.secondary_yaxis('right',
                                functions=(functools.partial(values_counts_to_mag,const=const),
                                            functools.partial(values_mag_to_counts,const=const)))
    
    counts_ticks = ax.get_yticks()
    mag_ticks = values_counts_to_mag(counts_ticks,const=const)
    axyrig.set_yticks(mag_ticks)
    axyrig.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axyrig.yaxis.set_minor_formatter(plt.NullFormatter())
    
    axyrig.set_ylabel(r'$\ln(\mu/\mathrm{mag/arcsec^2})$')
    #axyrig.tick_params(axis='y', which='major') 
    
    #axyrig.minorticks_on()
    
    if '-png' in sys.argv:
        img_format = 'png'
    else: 
        img_format = 'pdf'
    
    fig_path = f'{galaxy_folder_path}/{fig_name}'
    plt.savefig(f'{fig_path}.{img_format}', format=img_format, dpi=1000, bbox_inches='tight')   
    plt.close()