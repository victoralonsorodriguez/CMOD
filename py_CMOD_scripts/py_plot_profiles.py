
import pandas as pd
import numpy as np
import functools
import sys

from astropy.cosmology import FlatLambdaCDM


import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from cmod.photometry import values_counts_to_mag,values_mag_to_counts
from cmod.utils import rad_to_deg, deg_to_rad, ell_to_axrat,axrat_to_ell
from py_convert_functions import px_to_kpc,kpc_to_px
                                  



def plot_profiles(csv_path_list='.',
                  fig_name='fig_name',
                  cons=None,
                  final_plot=False,
                  output_path='.'):
    
    
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(0.005987)

    
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
    
    
    csv_values_dict = {}
    df_len_list = []
    for csv_label in csv_path_list:
        df = pd.read_csv(csv_label[0])
        csv_values_dict[csv_label[1]] = df
        df_len_list.append(len(df))

    max_data_len = min(df_len_list)

    profile_to_plot = ['ellipticity','pa','intens']
    profile_axis_label = ['Ellipticity',
                          'Position Angle [deg]',
                          'Intensity [counts]']
    
    if final_plot == True:
        
        profile_to_plot = ['intens']
        profile_axis_label = ['Intensity [counts]']
        
    # Creating some figures
    plot_rows = 1
    plot_cols = len(profile_to_plot)
    
    fig = plt.figure(figsize=(6*plot_cols, 5*plot_rows))
    plt.subplots_adjust(hspace=0.5, wspace=0.5)

    colors = ['gold','deepskyblue','lime','deeppink']
    markers = ['o','v','*','X']
    
    # Every iteration of this loop is a new plot
    for prof_pos,prof in enumerate(profile_to_plot):
        color_index = 0
        y_max = -10**20
        y_min = 10**20 
        ax = plt.subplot(plot_rows, plot_cols, prof_pos+1)
        
        for plot_label in csv_values_dict.keys():

            if plot_label == 'func_df':
                total_int_label_list = []
                for col in csv_values_dict[plot_label].columns:
                        if col != 'sma' and col != 'Total_int':
                            total_int_label_list.append(col)  
                result = ' + '.join(f'{label}' for label in total_int_label_list)

                for col in csv_values_dict[plot_label].columns:
        
                    if col != 'sma':
                        
                        x = csv_values_dict[plot_label]['sma'][:max_data_len]
                        y = csv_values_dict[plot_label][f'{col}'][:max_data_len]
                        
                        if max(y) > y_max:
                            y_max = max(y)
                        if min(y) < y_min:
                            y_min = min(y)

                        if col == 'Total_int':
                            
                            plt.plot(x,y,label=f'{result}',
                                    linewidth=1,color='black',
                                    linestyle='-.',zorder=4)
                            
                        else:
                            plt.scatter(x,y,label=col,
                                        marker=markers[color_index],s=10,
                                        linewidth=0.15,edgecolor='black',
                                        color=colors[color_index],zorder=1)
                            color_index += 1

            else:
                
                x = csv_values_dict[plot_label]['sma'][:max_data_len]
                y = csv_values_dict[plot_label][f'{prof}'][:max_data_len]
                
                # Isohpohes pos angle is in rad so can be changed to deg
                if profile_axis_label[prof_pos] == 'Position Angle [deg]':
                    angle_deg = rad_to_deg(y)   
                    y = angle_deg      
                    y_error = csv_values_dict[plot_label][f'{prof}_err'] * 180 / np.pi
                
                if max(y) > y_max:
                        y_max = max(y)
                if min(y) < y_min:
                        y_min = min(y)
                
                z_order = 2
                marker_size = 11
                if plot_label == 'Model':
                    z_order = 3
                    marker_size = 14
                plt.scatter(x,y,label=plot_label,
                            marker=markers[color_index],s=marker_size,
                            linewidth=0.15,edgecolor='black',
                            color=colors[color_index],zorder=z_order)
                color_index += 1

        # Legend position
        plt.legend(loc='best',prop={'size': 10})

        # Customizing the plots       
        # X bottom axis is common
        ax.set_xlabel(r'$Semimajor\ Axis\ Length\ [\mathrm{pix}]$')
        
        ax.set_xticks(np.linspace(np.min(x), np.max(x), 7))
        ax.set_xmargin(0.1)
        ax.minorticks_on()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        
        ax.grid(True,alpha=0.2)
        
        # TOP X axis is common
        axxtop = ax.secondary_xaxis('top',
                                    functions=(functools.partial(px_to_kpc,inst=cons[0]),
                                               functools.partial(kpc_to_px,inst=cons[0])))
        
        px_ticks = ax.get_xticks()
        arcsec_ticks = px_to_kpc(px_ticks,inst=cons[0])
        axxtop.set_xticks(arcsec_ticks)
        
        axxtop.minorticks_on()
        
        axxtop.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        
        axxtop.set_xlabel(r'$Semimajor\ Axis\ Length\ [\mathrm{kpc}]$',labelpad=8)
        axxtop.tick_params(axis='x', which='major')
        
        # Y axes are different for each plot
        # For intensisty
        if prof == 'intens':
        
            # Y left axis
            ax.set_ylabel(r'$Intensity\ [\mathrm{counts}]$')
            ax.set_yticks(np.linspace(np.min(y_min), np.max(y_max), 7))
            ax.minorticks_on()
            
            ax.set_yscale('log')
            #ax.set_ylim(bottom=0.5)
            
            # Y right axis
            axyrig = ax.secondary_yaxis('right',
                                        functions=(functools.partial(values_counts_to_mag,const=cons),
                                                   functools.partial(values_mag_to_counts,const=cons)))
            
            counts_ticks = ax.get_yticks()
            mag_ticks = values_counts_to_mag(counts_ticks,const=cons)
            axyrig.set_yticks(mag_ticks)
            axyrig.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axyrig.yaxis.set_minor_formatter(plt.NullFormatter())
            
            axyrig.minorticks_on()
            
            axyrig.set_ylabel(r'$\mu\ [\mathrm{mag/arcsec^2}]$')
            axyrig.tick_params(axis='y', which='major')
            
        elif prof == 'ellipticity':
            
            # Y left axis
            ax.set_ylabel(r'$Ellipticity$')
            ax.set_yticks(np.linspace(np.min(y_min), np.max(y_max), 7))
            ax.minorticks_on()
            
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            
            # Y right axis
            axyrig = ax.secondary_yaxis('right',functions=(ell_to_axrat,axrat_to_ell))
            
            ell_ticks = ax.get_yticks()
            axrat_ticks = ell_to_axrat(ell_ticks)
            axyrig.set_yticks(axrat_ticks)
            axyrig.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axyrig.yaxis.set_minor_formatter(plt.NullFormatter())
            
            axyrig.minorticks_on()
            
            axyrig.set_ylabel(r'$Axis\ ratio$')
            axyrig.tick_params(axis='y', which='major')
            
            
        elif prof == 'pa':
        
            # Y left axis
            ax.set_ylabel(r'$Position\ Angle\ [\mathrm{deg}]$')
            ax.set_yticks(np.linspace(np.min(y_min), np.max(y_max), 7))
            ax.minorticks_on()
            
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            
            # Y right axis
            axyrig = ax.secondary_yaxis('right',functions=(deg_to_rad,rad_to_deg))
            
            deg_ticks = ax.get_yticks()
            rad_ticks = deg_to_rad(deg_ticks)
            axyrig.set_yticks(list(rad_ticks))
            axyrig.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axyrig.yaxis.set_minor_formatter(plt.NullFormatter())
            
            axyrig.minorticks_on()
            
            axyrig.set_ylabel(r'$Position\ Angle\ [\mathrm{rad}]$')
            axyrig.tick_params(axis='y', which='major')
            
            
    if '-png' in sys.argv:
        img_format = 'png'
    else: 
        img_format = 'pdf'

    fig_name_final = f'{fig_name}_image_profiles.{img_format}'
    fig_path = f'{output_path}/{fig_name_final}'
    plt.savefig(f'{fig_path}', format=img_format, dpi=1000, bbox_inches='tight')    
    plt.close()    