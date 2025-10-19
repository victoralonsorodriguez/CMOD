
import numpy as np
import pandas as pd
pd.set_option('mode.chained_assignment', None)
import sys
import pdb
import functools
from itertools import groupby
import os
import copy
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FixedLocator
import matplotlib.font_manager as font_manager

from astropy.cosmology import FlatLambdaCDM

from .cosmology import px_to_kpc,kpc_to_px, RtW, WtR
from .io import open_fits, argparse_values
from .photometry import filter_wl_dict,values_counts_to_mag, values_mag_to_counts
from .utils import (rad_to_deg, deg_to_rad, ell_to_axrat, axrat_to_ell, 
                    create_folder, version_file_last, Cronometro)


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



# Enable LaTeX rendering
# Define the LaTeX preamble with siunitx
plt.rcParams['text.latex.preamble'] = r'''
            \usepackage{siunitx}
            \usepackage{bm}
            \sisetup{
            detect-family,
            separate-uncertainty=true,
            output-decimal-marker={.},
            exponent-product=\cdot,
            inter-unit-product=\cdot,
            }
            \DeclareSIUnit{\cts}{cts}
            '''
            
plt.rc('text', usetex=True)
plt.rc('font', family='serif')  # Use a serif font for LaTeX rendering




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
    
    
    
    


# Position of the minor ticks in top axis
def minor_ticks_pos(z_labels,filter_value):
    tick_pos_list = []
    
    for i in range(len(z_labels)-1):
        tick_pos_z = np.linspace(z_labels[i],z_labels[i+1],5)
        for i in range(len(tick_pos_z)-1):
            tick_pos_list.append(RtW(tick_pos_z[i],filter_value))

    return np.array(tick_pos_list)

def list_of_tuple(lst):
    if isinstance(lst, str):  
        lst = [lst]
    return [tuple(lst)]


# Creating a general function to plot
def plot_gen(subplots_conf = False,
             subplots_num = (1,1),
             subplots_positions = [(0,0)],
             subplots_titles = None,
             subplots_split_groups = 1,
             subplots_save_ind = False,
             x_data=[],
             y_data=[],
             label_list=None,
             legend_plot = True,
             legend_pos = 'best',
             legend_font_size  = 16,
             leg_col=1,
             fig_fonsize = 24,
             fig_title = None,
             fig_name='figure_withnoname',
             fig_save_path = '.',
             fig_format = 'pdf',
             aspect_ratio=[1,1],
             back_color = None,
             axis_font_size = 24,
             axix_grid = True,
             x_axis_label = 'X label',
             x_con_factor=1,
             x_log_scale = False,
             x_invert=False,
             x_lim=None,
             x_ticks_lim = None,
             x_ticks_num = 5,
             x_ticks_dec = 2,
             x_secaxix_ticks = False,
             x_axis_sides = 'bottom',
             x_top_axis_dec = 1,
             x_diff_fuct = None,
             y_axis_label = 'Y label',
             y_con_factor=1,
             y_log_scale = False,
             y_invert = False,
             y_lim = None,
             y_ticks_lim = None,
             y_ticks_num = 5,
             y_ticks_dec = 2,
             y_axis_color = False,
             y_axis_diff = False,
             y_axis_sides = 'left',
             y_diff_fuct_val = None,
             legend_pos_diff = 'lower right',
             y_secaxix_ticks = False,
             line_color = None,
             line_alpha = None,
             line_style = None,
             width_style = None,
             mark_style = None,
             mark_size = None,
             zorder = None,
             scatter_plot = False,
             style_plot = [], # 'sct', 's-l' or 'line'
             point_lab=[False,[],[[None,None]],None,[],None], # Yes/No,labels,text_pos,fontsize,selected_data,color
             tex_anno=[False,[],None,None], # Yes/No,labels,text_pos,fontsize,selected_data
             guide_lines = [False,(None,None)],
             fill_area=[False,[]],
             plot_show = False):
    
    plt.rc('font', size=fig_fonsize)  # Adjust size to your preference
    
    
    # Determinando el tamaño de las figuras
    
    # Cuánto miden de alto
    if aspect_ratio != 'auto':
        plot_cols = 1
        size_factor_cols = 8*aspect_ratio[0]
        fig_size_cols = size_factor_cols * plot_cols
    
        # Cuánto miden de ancho
        plot_rows = 1
        size_factor_rows = 8*aspect_ratio[1]
        fig_size_rows = size_factor_rows * plot_rows

    # Tamaño de ls figuras
    fig_size_global = (fig_size_rows,fig_size_cols)
    
    # Leyenda
    if label_list == None:
        label_list = ['']*len(y_data)
    
    # Definiendo algunos colores
    if line_color is None:
        line_color_list = [('red','dodgerblue','lime','blueviolet','darkorange','black','deeppink','gold')]
        line_color = line_color_list*(len(x_data)//len(line_color_list) + 1)

    # Definiendo alpha:
    if line_alpha is None:
        line_alpha = [1]*len(x_data)
        
    # zorder
    if zorder is None:
        zorder = [3]*len(x_data)
        if y_axis_sides == 'both':
            zorder = zorder*len(x_data[0])

    # Estilos de linea
    if line_style is None:
        line_style = ['-', '--', '-.']*(len(x_data)//3 + 1)
    
    # Grosor de linea
    if width_style is None:
        width_style = [1]*len(x_data)

    # Estilos de marcadores de scatter
    if mark_style is None:
        mark_style_list = ['o', '<','*','>', 'v']
        mark_style = mark_style_list*(len(x_data)//len(mark_style_list) + 1)
    
    if mark_size is None:
        mark_size = [20,20,50,20,20]
        mark_size = mark_size*(len(x_data)//len(mark_size) + 1)
        
    if len(style_plot) == 0:
        if scatter_plot == False:
            style_plot = ['line']*len(x_data)
        else:
            style_plot = ['sct']*len(x_data)
            
    elif len(style_plot)>0 and len(style_plot)<len(x_data):
        if scatter_plot == False:
            style_plot = style_plot*(len(x_data)//len(style_plot) + 1)
        else:
            style_plot = ['sct']*(len(x_data)//len(style_plot) + 1)
    

    if y_axis_sides != 'both':

        x_axis_label = list_of_tuple(x_axis_label)    
        y_axis_label = list_of_tuple(y_axis_label)    
        label_list = list_of_tuple(label_list)    
        line_color = list_of_tuple(line_color)    
        line_alpha = list_of_tuple(line_alpha)    
        zorder = list_of_tuple(zorder)    
        line_style = list_of_tuple(line_style)    
        width_style = list_of_tuple(width_style)   
        mark_style = list_of_tuple(mark_style) 
        mark_size = list_of_tuple(mark_size) 
        style_plot = list_of_tuple(style_plot)  
        #y_invert = list_of_tuple(y_invert)
    
    # Informando en la terminal sobre lo que se está ploteando
    print(f'\nPlotting: {fig_name.split("/")[-1]}\n')
    
    # Creando la figura como un solo subplot con ejes 'ax'
    #fig = plt.figure(figsize=fig_size_global)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    sub_row = subplots_num[0]
    sub_col = subplots_num[1]
    
    # Para cambiar las propiedades de los plots cada 
    # cierto número de filas
    subplots_split_row_num = sub_row // subplots_split_groups
    
    fig, axs = plt.subplots(nrows=sub_row, ncols=sub_col,
                    figsize=(sub_col*fig_size_rows, sub_row*fig_size_cols))
    
    # Remove empty subplots
    for row in range(sub_row):
        for col in range(sub_col):
            if (row, col) not in [pos for pos in subplots_positions]:
                if sub_col > 1:
                    fig.delaxes(axs[row, col])
                else:
                    fig.delaxes(axs[row])
                
    sub_plot_cost = True
    
    for data_pos in range(len(x_data)):
        
        if subplots_conf == True:
            row = subplots_positions[data_pos][0] 

            # Dividing the plots into gropus with different properties
            data_group = row // subplots_split_row_num 
            
            col = subplots_positions[data_pos][1]
        
        else:
            data_group = 0
    
        if subplots_save_ind == True:
            fig, ax = plt.subplots(nrows=1, ncols=1,
                    figsize=(1*fig_size_rows, 1*fig_size_cols))
            
        else:
            # If there is just one plot
            if sub_col == 1 and sub_row > 1:
                ax = axs[row]
            
            # For more than one plot
            elif sub_col > 1:
                ax = axs[row][col]
                
            elif sub_col == 1 and sub_row == 1:
                ax = axs
            
        
        # Different index for data and plot atributes
        if y_axis_sides == 'both':
            prop_pos = data_pos % 2
        else:
            prop_pos = data_pos
            
        prop_pos_l = 0
        prop_pos_r = 1
            
        
        # Si se indica un título se pone. Por defecto desactivado
        if '--no-title' not in sys.argv:
            if fig_title != None:
                fig.suptitle(f'{fig_title}',fontsize=axis_font_size+10,y=0.995)
        
        if subplots_titles != None:
            ax.set_title(f'\\textbf{{{subplots_titles[data_pos]}}}',fontsize=axis_font_size+1)
            
        
        x_data_bottom = np.array(x_data[data_pos])/x_con_factor
        y_data_left = np.array(y_data[data_pos])/y_con_factor
        
        if y_axis_sides == 'both':
            y_data_left = np.array(y_data[data_pos][0])/y_con_factor
            y_data_right = np.array(y_data[data_pos][1])/y_con_factor
        
        # Representando los datos
        if y_axis_sides == 'both' or y_axis_sides == 'right':
            axyrig = ax.twinx()
        
        if style_plot[data_group][prop_pos] == 'sct':

            if y_axis_sides == 'left' or y_axis_sides == 'both':
                if y_axis_sides == 'both':
                    prop_pos = prop_pos_l
                
                ax.scatter(x_data_bottom,
                        y_data_left,
                        label=f'{label_list[data_group][prop_pos]}',
                        linewidth=0.15,
                        edgecolor='black',
                        color=line_color[data_group][prop_pos],
                        s=mark_size[data_group][prop_pos],
                        marker=mark_style[data_group][prop_pos],
                        zorder=zorder[data_group][prop_pos],
                        alpha=line_alpha[data_group][prop_pos])  
                
            if y_axis_sides == 'right' or y_axis_sides == 'both':  
                if y_axis_sides == 'both':
                    prop_pos = prop_pos_r
                axyrig.scatter(x_data_bottom,
                        y_data_right,
                        label=f'{label_list[data_group][prop_pos]}',
                        linewidth=0.15,
                        edgecolor='black',
                        color=line_color[data_group][prop_pos],
                        s=mark_size[data_group][prop_pos],
                        marker=mark_style[data_group][prop_pos],
                        zorder=zorder[data_group][prop_pos],
                        alpha=line_alpha[data_group][prop_pos]) 
                        
            
        elif style_plot[data_group][data_pos] == 's-l':
            if zorder[data_group][data_pos] == 0:
                zorder[data_group][data_pos] = 1
                
            if y_axis_sides == 'left' or y_axis_sides == 'both':    
                ax.scatter(x_data_bottom,
                        y_data_left,
                        label=f'{label_list[data_group][data_pos]}',
                        linewidth=0.15,
                        edgecolor='black',
                        color=line_color[data_group][prop_pos],
                        s=mark_size[data_group][prop_pos],
                        marker=mark_style[data_group][prop_pos],
                        zorder=zorder[data_group][prop_pos])

                ax.plot(x_data_bottom,
                        y_data_left,
                        linewidth=0.5,
                        color=line_color[data_group][prop_pos],
                        linestyle=line_style[data_group][prop_pos],
                        zorder=zorder[data_group][prop_pos]-1)
                
            if y_axis_sides == 'right' or y_axis_sides == 'both':
                axyrig.scatter(x_data_bottom,
                        y_data_right,
                        label=f'{label_list[data_group][data_pos]}',
                        linewidth=0.15,
                        edgecolor='black',
                        color=line_color[data_group][data_pos+1],
                        s=mark_size[data_group][data_pos+1],
                        marker=mark_style[data_group][data_pos+1],
                        zorder=zorder[data_group][data_pos+1])

                axyrig.plot(x_data_bottom,
                        y_data_right,
                        linewidth=0.5,
                        color=line_color[data_group][data_pos+1],
                        linestyle=line_style[data_group][data_pos+1],
                        zorder=zorder[data_group][data_pos+1]-1)
                
            
        elif style_plot[data_group][data_pos] == 'line':
            if y_axis_sides == 'left' or y_axis_sides == 'both':  
                ax.plot(x_data_bottom,
                        y_data_left,
                        label=f'{label_list[data_group][data_pos]}',
                        color=line_color[data_group][prop_pos],
                        linestyle=line_style[data_group][prop_pos],
                        linewidth=width_style[data_group][prop_pos],
                        zorder=zorder[data_group][prop_pos],
                        alpha=line_alpha[data_group][prop_pos])
            if y_axis_sides == 'right' or y_axis_sides == 'both':  
                axyrig.plot(x_data_bottom,
                        y_data_left,
                        label=f'{label_list[data_group][data_pos+1]}',
                        color=line_color[data_group][data_pos+1],
                        linestyle=line_style[data_group][data_pos+1],
                        linewidth=width_style[data_group][data_pos+1],
                        zorder=zorder[data_group][data_pos+1],
                        alpha=line_alpha[data_group][data_pos+1])



        # Para rellenar los areas
        if fill_area[0] == True:
            fill_max = max(max(row) for row in y_data)
            if len(fill_area[1]) > 0:
                if data_pos+1 in fill_area[1]:
                    ax.fill_between(x_data[data_pos],y_data[data_pos],fill_max,
                                    color='black', alpha=fill_area[2],
                                    lw = 0)
                    
                if data_pos + 1 == 9:
                    ax.fill_between(x_data[data_pos],y_data[data_pos],fill_max,
                                    color='black', alpha=fill_area[2],
                                    lw = 0)
                    ax.fill_between(x_data[data_pos-1],y_data[data_pos-1],fill_max,
                                    color='black', alpha=fill_area[2],
                                    lw = 0)
            else:
                ax.fill_between(x_data[data_pos],y_data[data_pos],fill_max,
                                    color='black', alpha=fill_area[2],
                                    lw = 0)
        

        
        if point_lab[0] == True:
            # Text position in a tuple at 3rd pos
            if len(point_lab)<3 or (len(point_lab)>=3 and len(point_lab[2])==0):
                point_lab_x = 0
                point_lab_y = 0
            elif len(point_lab)>=3:
                if len(point_lab[2]) == 0:
                    point_lab_x = 0
                    point_lab_y = 0   
                elif len(point_lab[2]) > 0:
                    if len(point_lab[2]) == 1:
                        point_lab_x = point_lab[2][0][0]
                        point_lab_y = point_lab[2][0][1]  
                    if len(point_lab[2]) > 1 and data_pos in np.arange(len(point_lab[2])):
                        if style_plot[data_pos] == 'line':
                            point_lab_x = point_lab[2][data_pos][0]
                            point_lab_y = point_lab[2][data_pos][1] 
                
            if len(point_lab)<4 or point_lab[3] == None:
                point_lab_size = 15
            elif len(point_lab)>=4 and point_lab[3] != None:
                point_lab_size = point_lab[3]
            
            if len(point_lab)<6 or point_lab[5] == None:
                lab_color = 'black'
            elif len(point_lab)>=6 and point_lab[5] != None:
                lab_color = point_lab[5]
            
            if style_plot[data_group][data_pos] == 'sct':
                if len(point_lab)<5 or (len(point_lab)>=5 and len(point_lab[4])==0):
                    for tex_pos in range(len(x_data[0])):
                        point_lab_x = point_lab[2][tex_pos][0]
                        point_lab_y = point_lab[2][tex_pos][1]    
                        ax.annotate(point_lab[1][tex_pos], (x_data[data_pos][tex_pos], y_data[data_pos][tex_pos]),
                                    fontsize=point_lab_size,xytext=(point_lab_x,point_lab_y),textcoords='offset pixels',
                                    color=lab_color)
                        
                elif len(point_lab)>=5 and len(point_lab[4])>0:
                    if data_pos+1 in point_lab[4]:
                        for tex_pos in range(len(x_data[0])):
                            point_lab_x = point_lab[2][tex_pos][0]
                            point_lab_y = point_lab[2][tex_pos][1]  
                            ax.annotate(point_lab[1][tex_pos], (x_data[data_pos][tex_pos], y_data[data_pos][tex_pos]),
                                    fontsize=point_lab_size,xytext=(point_lab_x,point_lab_y),textcoords='offset pixels',
                                    color=lab_color)
                    else:
                        pass
                    
            if style_plot[data_group][data_pos] == 'line':
                if len(point_lab)<5 or (len(point_lab)>=5 and len(point_lab[4])==0):
                        ax.annotate(point_lab[1][data_pos], (x_data[data_pos][0], y_data[data_pos][0]),
                                    fontsize=point_lab_size,xytext=(point_lab_x,point_lab_y),textcoords='offset pixels',
                                    color=lab_color)
                        
                elif len(point_lab)>=5 and len(point_lab[4])>0:
                    if data_pos+1 in point_lab[4]:
                            ax.annotate(point_lab[1][data_pos], (x_data[data_pos][0], y_data[data_pos][0]),
                                    fontsize=point_lab_size,xytext=(point_lab_x,point_lab_y),textcoords='offset pixels',
                                    color=lab_color)
                    else:
                        pass  
                        
        if tex_anno[0] == True:
            if len(tex_anno)<3 or tex_anno[2] == None:
                    point_tex_size = 10
            elif len(tex_anno)>=3 or tex_anno[2] != None:
                    point_tex_size = tex_anno[2]
                    
            if len(tex_anno)<4 or tex_anno[3] == None:
                    tex_col = 'black'
            elif len(tex_anno)>=4 or tex_anno[3] != None:
                    tex_col = tex_anno[3]
                
            for pos in range(len(tex_anno[1])):
                text = tex_anno[1][pos][0]            
                text_x_pos,text_y_pos = tex_anno[1][pos][1]
                if len(tex_anno[1][pos])>=3:
                    rot_ang=tex_anno[1][pos][2]
                else:
                    rot_ang = 0
                plt.annotate(text, (text_x_pos,text_y_pos),
                            fontsize=point_tex_size,xytext=(0,0),textcoords='offset pixels',
                            color=tex_col,rotation=rot_ang)

        

        if back_color != None:
            ax.set_facecolor(back_color)    
                            
        # Para añadir lineas verticales de referencia
        # Por defecto desactivadas
        if guide_lines[0]==True:
            
            if guide_lines[1][0] != None and len(guide_lines[1][0])>0:
                for line_pos in guide_lines[1][0]:
                    ax.axvline(line_pos, color="black",linestyle="dashdot", linewidth=1, 
                        zorder=0, alpha=0.9)

            if guide_lines[1][1] != None and len(guide_lines[1][1])>0:
                for line_pos in guide_lines[1][0]:
                    ax.axhline(line_pos, color="black",linestyle="dashdot", 
                               linewidth=1, zorder=0, alpha=0.9)
                    
    
            
        ######################################################################    
        ###-----------------------AXIS CUSTOMIZATION-----------------------###
        ######################################################################
        
        
        
        ###------------------X BOTTOM AXIS------------------###
        #######################################################
        
        ax.set_xlabel(f'{x_axis_label[data_group][0]}',fontsize=axis_font_size)

        # Si se limita el eje X
        if x_lim != None and x_lim != 'plt' and len(x_lim)>0:
                if len(x_lim) == 2:
                    ax.set_xlim(left=min(x_lim), right=max(x_lim))
                elif len(x_lim) == 1 and isinstance(x_lim[0], (int, float)):
                    ax.set_xlim(left=x_lim[0])
                elif len(x_lim) == 1 and isinstance(x_lim[1], (int, float)):
                    ax.set_xlim(right=x_lim[1])
        
        
        
        # Si se dan los limites de ticks se ponen sino se calculan
        if x_ticks_lim != 'plt':
            if x_ticks_lim != None:
                ax.set_xticks(np.linspace(min(x_ticks_lim), max(x_ticks_lim), x_ticks_num))
            elif x_lim != None and x_lim != 'plt' and x_lim != None:
                ax.set_xticks(np.linspace(min(x_lim), max(x_lim), x_ticks_num))
            else:
                x_min = min(min(row) for row in x_data)
                x_max = max(max(row) for row in x_data)
                x_ticks_list = (np.linspace(x_min/x_con_factor, x_max/x_con_factor, x_ticks_num)).tolist()
                ax.set_xticks(x_ticks_list)

        ax.set_xmargin(0.05)
        ax.minorticks_on()
        
        # Numero de decimales de las etiquetas de los ticks
        ax.xaxis.set_major_formatter(FormatStrFormatter(f'%.{x_ticks_dec}f'))
        
        # Escala logarítmica si se indica
        if x_log_scale == True:
            ax.set_xscale('log')
            
        # Invertir el eje si se indica
        if x_invert == True:
            ax.invert_xaxis()

        ###------------------X TOP AXIS------------------###
        ####################################################
        
        if x_axis_sides == 'both':

            if x_diff_fuct != None:
                if x_diff_fuct[0] == 'RtW':
                
                    # Transforming redshift into wavelength
                    filter_wl_value = x_diff_fuct[1][data_pos]
                    axxtop = ax.secondary_xaxis('top', 
                            functions=(functools.partial(RtW, filter_wl=filter_wl_value),
                                        functools.partial(WtR, filter_wl=filter_wl_value)))
                    
                    axxbot_ticks = ax.get_xticks()
                    wl_ticks = RtW(axxbot_ticks,filter_wl_value)
                    axxtop.set_xticks(wl_ticks)

                    ticks_pos = minor_ticks_pos(axxbot_ticks,filter_wl_value)
                    axxtop.xaxis.set_minor_locator(FixedLocator(ticks_pos))
                    axxtop.set_ticks(ticks_pos, minor=True)
                    
                    axxtop.xaxis.set_major_formatter(FormatStrFormatter(f'%.{x_top_axis_dec}f'))
        
            axxtop.set_xlabel(f'{x_axis_label[data_group][1]}',fontsize=axis_font_size,labelpad=12)
            
        elif x_axis_sides == 'bottom' and sub_plot_cost == True:
            # Creamos el eje
            axxtop = ax.secondary_xaxis('top')
            
            # Obtenemos los ticks inferiores
            axxbot_ticks = ax.get_xticks()
            
            # Para poner los labels sobre los ticks, por defecto se ponen
            if x_secaxix_ticks == False:
                axxtop.xaxis.set_major_formatter(FormatStrFormatter(''))
        
            axxtop.xaxis.set_major_locator(FixedLocator(axxbot_ticks))
            
            minor_ticks = ax.xaxis.get_minorticklocs()
            axxtop.xaxis.set_minor_locator(FixedLocator(minor_ticks))
            
            # Numero de decimales en los ticks
            axxtop.xaxis.set_major_formatter(FormatStrFormatter(f'%.{x_ticks_dec}f'))
    
        
        # Escala logarítmica si se indica
        if x_log_scale == True:
            axxtop.set_xscale('log')
            
        # Invertir el eje si se indica
        if x_invert==True:
            axxtop.invert_xaxis()



        ###------------------Y LEFT AXIS------------------###
        ####################################################
        
        if y_axis_color == True:
            y_left_axis_color = line_color[data_group][prop_pos_l]
        else:
            y_left_axis_color = 'black'
        
        if y_axis_sides == 'both':
            ax.set_ylabel(y_axis_label[data_group][prop_pos_l],fontsize=axis_font_size,
                        labelpad=10, color=y_left_axis_color)
        else:
            ax.set_ylabel(y_axis_label[data_group][prop_pos_l],fontsize=axis_font_size,labelpad=10)
        
        # Si se limita el eje Y
        if y_lim != None and y_lim != 'plt' and len(y_lim)>0:
                if len(y_lim) == 2:
                    ax.set_ylim(bottom=min(y_lim), top=max(y_lim))
                elif len(y_lim) == 1 and isinstance(y_lim[0], (int, float)):
                    ax.set_ylim(bottom=y_lim[0])
                elif len(y_lim) == 1 and isinstance(y_lim[1], (int, float)):
                    ax.set_ylim(top=y_lim[1])
                    
        # Datos para plotear
        y_axis_l_data = np.copy(y_data_left)
        if y_axis_sides == 'both':
            y_axis_r_data = np.copy(y_data_right)   

        # Si tenemos eje secundario diferente
        #if y_axis_diff == True:
            
            # Eliminamos los valores del eje derecho
            #del y_axis_l_data[y_axis_diff_val_pos-1]
        
        # Si se dan los limites de ticks se ponen sino se calculan
        if y_ticks_lim != 'plt' and sub_plot_cost == True:
            if y_ticks_lim != None:
                ax.set_yticks(np.linspace(min(y_ticks_lim), max(y_ticks_lim), y_ticks_num))
            elif y_lim != None and y_lim != 'plt' and y_lim != None:
                ax.set_yticks(np.linspace(min(y_lim), max(y_lim), y_ticks_num))
            else:
                #y_min = min(min(row) for row in y_axis_l_data)
                #y_max = max(max(row) for row in y_axis_l_data)
                y_min = np.nanmin(y_axis_l_data)
                y_max = np.nanmax(y_axis_l_data)
                y_ticks_l_list = (np.linspace(y_min/y_con_factor, y_max/y_con_factor, y_ticks_num)).tolist()
                ax.set_yticks(y_ticks_l_list)
                    
        #ax.set_ymargin(0.1)
        ax.yaxis.set_major_formatter(FormatStrFormatter(f'%.{y_ticks_dec}f'))
        ax.yaxis.set_minor_formatter(plt.NullFormatter())   
        #ax.minorticks_on()

        # Escala logarítmica si se indica
        if y_log_scale == True:
            ax.set_yscale('log')
        
        # Invertir el eje si se indica
        if y_invert==True:
            ax.invert_yaxis()

        # Color para los ejes
        if y_axis_color == True:
            ax.tick_params(axis='y', labelcolor=y_left_axis_color)

        ###------------------Y RIGHT AXIS------------------###
        ####################################################
        
        # Si no tenemos eje secundario diferente creamos uno
        # que es igual que el izquierdo
        if y_axis_color  == True:
            y_right_axis_color = line_color[data_group][prop_pos_r]
        else:
            y_right_axis_color = 'black'
        
        if y_axis_diff == False:
            # Si no tenemos eje secundario seguimos como siempre
            axyrig = ax.secondary_yaxis('right')

        else:
            
            # Nombre del eje secundario diferente
            axyrig.set_ylabel(f'{y_axis_label[data_group][prop_pos_r]}',fontsize=axis_font_size,labelpad=10,
                            color=y_right_axis_color)
        
        # Obteniendo los ticks
        # Si no tenemos eje secundario diferente, como siempre
        if y_axis_diff == False:
            axylef_ticks = ax.get_yticks()
            axyrig.set_yticks(axylef_ticks)
        
        else:
            # Si tenemos un eje secundario diferente
            # volvemos a calcular los ticks
            y_min = np.nanmin(y_data_right)
            y_max = np.nanmax(y_data_right)
            y_ticks_r_list = (np.linspace(y_min/y_con_factor, y_max/y_con_factor, y_ticks_num)).tolist()
            axyrig.set_yticks(y_ticks_r_list)
            axyrig.yaxis.set_major_locator(FixedLocator(y_ticks_r_list))
            #axyrig.set_ymargin(0.1)

        
        axyrig.yaxis.set_major_formatter(FormatStrFormatter(f'%.{y_ticks_dec}f'))
        
        if y_axis_diff == True:
            y_secaxix_ticks = True
        
        # Para no poner los labels sobre los ticks, por defecto se ponen
        if y_secaxix_ticks == False:
            axyrig.yaxis.set_major_formatter(FormatStrFormatter(''))
    
        axyrig.yaxis.set_minor_formatter(plt.NullFormatter())   
        axyrig.minorticks_on()
        axyrig.set_ymargin(0.05)
            
        axyrig.tick_params(axis='y', which='major')
        
        # Escala logarítmica si se indica
        if y_log_scale == True:
            axyrig.set_yscale('log')
        
        # Invertir el eje si se indica
        
        #if subplots_conf == True:
        #    y_invert = line_color[data_group][prop_pos_l]
            
        if y_invert==True:
            axyrig.invert_yaxis()
            
        # Color para los ejes
        if y_axis_color == True:
            axyrig.tick_params(axis='y', labelcolor=y_right_axis_color)

        
        ###------------------OTHER PROPERTIES------------------###
        ##########################################################
        
        # Generando un grid sutil
        if axix_grid == True:
            ax.grid(True,alpha=0.2)  

        # Propiedades de la leyenda
        if legend_plot == True:
            font = font_manager.FontProperties(size=legend_font_size)
            ax.legend(loc=legend_pos,ncol=leg_col,prop=font,framealpha=1)
            if y_axis_diff == True:
                axyrig.legend(loc=legend_pos_diff,ncol=leg_col,prop=font,framealpha=1)
                
        # Dando espacio si se pone el título
        fig.tight_layout(pad=1)
    

    
        ###------------SAVING THE PLOT INDIVIDUALLY-----------###
        #########################################################   
        
        if subplots_save_ind == True:

            if subplots_conf == True:
                fig_path = f'{fig_save_path}/{fig_name}_{subplots_titles[data_pos]}'
            else:
                fig_path = f'{fig_save_path}/{fig_name}'
            
            plt.savefig(f'{fig_path}.{fig_format}',format=fig_format, dpi=1000, bbox_inches='tight')   
            
            if plot_show == True:
                plt.show()
                
            plt.close()
        
        if subplots_conf == False:    
            sub_plot_cost = False             
    
    ###------------------SAVING THE PLOT------------------###
    #########################################################   
    if subplots_save_ind == False:

        # Guardando los archivos en pdf o png desde la terminal.
        # Por defecto en pdf, para guardar en png usar -png en terminal
        fig_path = f'{fig_save_path}/{fig_name}'
        
        plt.savefig(f'{fig_path}.{fig_format}',format=fig_format, dpi=1000, bbox_inches='tight')   
            
        if plot_show == True:
            plt.show()
        
        # Cerrando el plot
        plt.close()   
    
    

def plot_prop_conf(mag_pairs_plot):
    
    mag_plot_prop_dict = {}
    
    # Plots configuration
    for mag_pair in mag_pairs_plot:
        for mag in mag_pair:
            
            mag_plot_style = 'sct'
            mag_legend_label = ' '
            mag_alpha = 0.75
            mag_zorder = 3
            mag_Y_ax_inversion = False
            
            if mag == 'Ax_rat':
                mag_label = 'Axis Ratio'
                mag_color = 'limegreen'
                mag_marker = '<'
                mag_marker_size = 150
                
            elif mag == 'Eff_rad_arcsec':
                mag_label = 'Effective Radius [arcsec]'
                mag_color = 'dodgerblue'
                mag_marker = '*'
                mag_marker_size = 200
                
            elif mag == 'Eff_rad_kpc':
                mag_label = 'Effective Radius [kpc]'
                mag_color = 'dodgerblue'
                mag_marker = '*'
                mag_marker_size = 200
            
            elif mag == 'Ind':
                mag_label = 'Sersic Index'
                mag_color = 'red'
                mag_marker = '>'
                mag_marker_size = 150
                
            elif mag == 'Mag':
                mag_label = 'Total Magnitude'
                mag_color = 'darkorange'
                mag_marker = 'o'
                mag_marker_size = 100
                mag_Y_ax_inversion = True
            
            elif mag == 'I_e':
                mag_label = 'Effective Intensity [counts]'
                mag_color = 'darkorange'
                mag_marker = 'o'
                mag_marker_size = 100    
            
            elif mag == 'Pos_ang':
                mag_label = 'Position Angle [°]'
                mag_color = 'gold'
                mag_marker = 'P'
                mag_marker_size = 150
                
            elif mag == 'X_cent':
                mag_label = 'Center X Position'
                mag_color = 'darkorchid'
                mag_marker = 'D'
                mag_marker_size = 100
                
            elif mag == 'Y_cent':
                mag_label = 'Center Y Position'
                mag_color = 'mediumturquoise'
                mag_marker = 'X'
                mag_marker_size = 150
            
            mag_plot_prop_dict[mag] = {}
            mag_plot_prop_dict[mag]['plot_style'] = mag_plot_style
            mag_plot_prop_dict[mag]['label'] = mag_label
            mag_plot_prop_dict[mag]['leg_lab'] = mag_legend_label
            mag_plot_prop_dict[mag]['color'] = mag_color
            mag_plot_prop_dict[mag]['marker'] = mag_marker
            mag_plot_prop_dict[mag]['marker_size'] = mag_marker_size
            mag_plot_prop_dict[mag]['alpha'] = mag_alpha
            mag_plot_prop_dict[mag]['zorder'] = mag_zorder
            mag_plot_prop_dict[mag]['Y_ax_inversion'] = mag_Y_ax_inversion
            
    return mag_plot_prop_dict



# Obtaining the key of a dictionary from a value
def get_key(val,dic):
    
    for key, value in dic.items():
        if val == value:
            return key


def filter_sort_dict(csv_folder_path,filter_wavelength_dict,filter_names_dict):
    
    # lists to store the filters names and the csv paths 
    filter_name_list = []
    file_path_list = []

    # Obtaining every csv file from the common created folder
    for file in sorted(os.listdir(csv_folder_path)):
                
        if '.csv' in file:
            
            # Removing the extension of the file
            csv_file_name = file.split('.')[0]

            # Obtaining the datacube name
            datacube_name = '_'.join(csv_file_name.split('_')[:-1])

            # Extracting the filter name
            filter_name = datacube_name.split('_')[-1]
            
            # Correcting somo filter names
            if filter_name[-1] == 'M':
                filter_name = filter_name.replace('M','W')
                
            if filter_name[:3] == 'HST':
                filter_name = filter_name.replace('W','H')
            
            # Adding filters to a list
            filter_name_list.append(filter_name)
            
            # Adding the absolute paths to a list
            file_path = os.path.join(csv_folder_path,file)
            file_path_list.append(file_path)    
        
    
    # Creating a sublist grouped by filter
    filter_name_group = [list(g) for _, g in groupby(filter_name_list, lambda k: k[-1])]

    # Correcting the JWT F###M filters name
    for i,group in enumerate(filter_name_group):
        for j,flt in enumerate(group):
            if any(str(flt[:-1]) in s for s in ('F140M','F162M','F182M','F210M')):
                filter_name_group[i][j] = f'{flt[:-1]}M'
                
    # Managing the HST filters name            
    for i,group in enumerate(filter_name_group):
        for j,flt in enumerate(group):
            if any(str(flt[:-1]) in s for s in ('HSTF300W', 'HSTF435W', 'HSTF450W',
                                                'HSTF475W', 'HSTF555W', 'HSTF569W', 
                                                'HSTF606W', 'HSTF791W', 'HSTF814W')):
                filter_name_group[i][j] = f'{flt[:-1]}W'            
                
    # Sorting the system group wavlenghts from min to max
    
    # To store the wavelength of each filter grouping by system
    filter_wavelength_group = [[] for i in range(len(filter_name_group))]
    
    # For the plot title
    filter_wavelength_group_title = [[] for i in range(len(filter_name_group))]

    # This loop sort the filters by alphabetical order and
    # the asociated wavelength in the same order
    for index_system,filter_system in enumerate(filter_name_group):
        for index_filtr,filtr in enumerate(filter_system):
            
            filter_wavelength_group[index_system].append(filter_wavelength_dict[filtr])
            filter_wavelength_group_title[index_system].append(filter_names_dict[filtr])

            filter_wavelength_group[index_system] = sorted(filter_wavelength_group[index_system])
            filter_wavelength_group_title[index_system] = sorted(
                                filter_wavelength_group_title[index_system])

    
    # For store de datacube filter names
    filter_name_group_sorted = [[] for i in range(len(filter_name_group))]
    
    # For store the correct order of filter oficial names for 
    # the plot title
    filter_name_group_title_sorted = [[] for i in range(len(filter_name_group))]

    # This loop sort the filters by increasing wavelength
    for index_wl_system,wl_system in enumerate(filter_wavelength_group):
        for index_wl_filtr,wl_filter in enumerate(wl_system):
            
            key = get_key(wl_filter,filter_wavelength_dict)
            key_title = filter_names_dict[key]

            filter_name_group_sorted[index_wl_system].append(key)
            filter_name_group_title_sorted[index_wl_system].append(key_title)

    df_filter_dict_large = {}
    df_filter_dict_short = {}   
    df_filter_dict_short_wl = {}
    df_filter_dict_large_to_short = {}
    
    # Creating an array plot for every filter system
    for filter_system_index, filter_system in enumerate(filter_name_group_sorted):
        
        # Differenciating the full name of the 
        # telescope and the acronym for the plot title
        if filter_system[0][:3] == 'Euc':
            filter_name_sys = 'Euclid'
            filter_name_title = 'Euclid'
        elif filter_system[0][0] == 'F':
            filter_name_sys = 'JamesWebb'
            filter_name_title = 'JWST'
        elif filter_system[0][:3] == 'HST':
            filter_name_sys = 'Hubble'
            filter_name_title = 'HST'
        elif filter_system[0][1:3] == 've':
            filter_name_sys = 'Johnson_UVRIJHK'
            filter_name_title = 'SDSS'
        
        df_index = 0
        
        # This dictionary contains as keys the filter SYSTEM names
        # Inside every SISTEM there is another dictionary with the
        # filter names as keys and the dataframes as values
        # We have a nested dictionary
        df_filter_dict_large[filter_name_title] = {}
        df_filter_dict_short[filter_name_title] = {}
        df_filter_dict_short_wl[filter_name_title] = {}
        
        for filter_name in filter_system:

            # Obtaining the csv path
            path_csv = [path for path in file_path_list if filter_name in path][0]
            
            # Storing the dataframes for each csv in a list
            df_filter_dict_large[filter_name_title][filter_name] = pd.read_csv(path_csv, sep=',')
            
            # For the short names
            filter_name_short = filter_name_group_title_sorted[filter_system_index][df_index]
            df_filter_dict_short[filter_name_title][filter_name_short] = pd.read_csv(path_csv, sep=',')
            
            # Pairing the short name with the wavelength
            filter_wl = filter_wavelength_dict[filter_name]
            df_filter_dict_short_wl[filter_name_title][filter_name_short] = filter_wl
            df_filter_dict_large_to_short[filter_name] = filter_name_short
            
            df_index = df_index + 1
            
    
    return df_filter_dict_short,df_filter_dict_short_wl,df_filter_dict_large_to_short



def plots_filter(system_list_mode,
                 filters_list_mode,
                 datacube_name_general,
                 analysis_programme,
                 galaxy_name,
                 plot_save_mode,
                 plot_max_columns,
                 plot_path,
                 filter_dict = {},
                 filters_to_plot = [],
                 mag_pairs_plot = [],
                 wl_dict = {},
                 plot_prop_dict = {}):     
    
    if system_list_mode == 0:
        filter_system_plot = filter_dict.keys()
    else:
        filter_system_plot = system_list_mode
    
    for mag_pair_pos,mag_pair in enumerate(mag_pairs_plot):
        
        # For storing the plot positions
        plots_positions = []
        
        # For storing the plots in the correct
        # order to plot
        plot_order = []
        
        # To avoiding duplicate filters
        filter_set = set()
        
        # For storing the filters names as 
        # titles of plots
        filter_suptitle = []
        
        # For storing the filters wavelenght
        # for convert from redshift to emitted wl
        filter_wl = []
        
        # Total columns of the plot, max is 6
        total_col  = 0
        
        # For store data in the correct order
        x_sub_data = []
        y_sub_data = []
        
        # Initizalizing the rows values
        # One row per filters system except if cols>6
        # for a given system
        rows = 0
        
        del_dict = copy.deepcopy(filter_dict)
        
        # For every filter system
        for filter_system in filter_system_plot:
            
            if filters_list_mode != 0:
                # Crear una lista con las claves a eliminar
                keys_a_eliminar = [key for key in del_dict[filter_system] if key not in filters_list_mode]
                # Eliminar las claves fuera del bucle
                for key in keys_a_eliminar:
                    del del_dict[filter_system][key]

            # Obtener las claves actualizadas después de la eliminación
            filters_to_plot = list(del_dict[filter_system].keys())
            if len(del_dict[filter_system].keys()) == 0:
                filters_to_plot = keys_a_eliminar            
            
            # Initializing column values
            cols = 0
            total_col_max = 0
            
            # For every filter in the filter system
            for filt in filters_to_plot:
                
                # If filter is in the list of filters to plot
                if filt in filters_to_plot and filt not in filter_set:
                    
                    plots_positions.append((rows,cols))
                    
                    # Avoiding to plot the same filter twice
                    filter_set.add(filt)
                        
                    # Name of each filter
                    filter_suptitle.append(filt)
                    
                    # Filter wavelength in the same order as the filters
                    filter_wl.append(wl_dict[filter_system][filt])
            
                    # Determining the order and the positions of the plots
                    plot_order.append(((filter_system,filt),(rows,cols)))
                    
                    # Obtainig the data for the plot
                    x_data = filter_dict[filter_system][filt]['z']
                    x_sub_data.append(x_data)
                    
                    # For the left Y axis
                    y_data_left = filter_dict[filter_system][filt][mag_pair[0]]
                    
                    # For the right Y axis
                    y_data_right = filter_dict[filter_system][filt][mag_pair[1]]
                    
                    # Tuple of y data values added to the y values list
                    y_sub_data.append((y_data_left,y_data_right))
                    
                    # Determining the number of columns for the plot
                    total_col_max += 1
                    if total_col_max > total_col:
                        total_col = total_col_max
                    
                    # Correction for systems with more than 
                    # max_col filters to plot by adding one extra row
                    if cols // plot_max_columns != 0 and cols == plot_max_columns:
                        rows += 1
                        cols = 0
                        plots_positions.pop()
                        plots_positions.append((rows,cols))
                    else:
                        cols += 1
                    
            # New row for each filter system
            rows += 1
        
        left_mag = mag_pair[0]
        left_mag_prop = plot_prop_dict[left_mag]
        
        right_mag = mag_pair[1]
        right_mag_prop = plot_prop_dict[right_mag]
        
        # Total rows for the plot sheet
        total_row = rows

        # Figure save name and figure title depending on:
        # If there is more than just one filter
        if len(filter_set) > 1 and plot_save_mode == 0:
            
            # If plotted mangitudes are the same 
            if left_mag_prop['label'] == right_mag_prop['label']:
                fig_name = f'{datacube_name_general}_{left_mag}'
                fig_title = (f'\\textbf{{{analysis_programme}}} output \\textbf{{{left_mag_prop["label"]}}} ' 
                            f'for galaxy \\textbf{{{galaxy_name}}}') 
                
            # If there are different magnitudes
            else:
                fig_name = f'{datacube_name_general}_{left_mag}_{right_mag}'
                fig_title = (f'\\textbf{{{analysis_programme}}} output \\textbf{{{left_mag_prop["label"]}}} and ' 
                        f'\\textbf{{{right_mag_prop["label"]}}} for galaxy \\textbf{{{galaxy_name}}}')
        
        # If there is only one filter
        else:
            
            # If plotted mangitudes are the same 
            if left_mag_prop["label"] == right_mag_prop["label"]:
                fig_name = f'{datacube_name_general}_{left_mag}'
                fig_title = (f'\\textbf{{{analysis_programme}}} output \\textbf{{{left_mag_prop["label"]}}}\n' 
                            f'for galaxy \\textbf{{{galaxy_name}}}') 
            else:
                fig_name = f'{datacube_name_general}_{left_mag}_{right_mag}'
                fig_title = (f'\\textbf{{{analysis_programme}}} output \\textbf{{{left_mag_prop["label"]}}} and \n' 
                        f'\\textbf{{{right_mag_prop["label"]}}} for galaxy \\textbf{{{galaxy_name}}}')
        
        if plot_save_mode == 0:
            subplots_save_ind = False
        else:
            subplots_save_ind = True

        # Plotting the filters with this configuration
        plot_gen(subplots_conf = True,
                subplots_num=(total_row,total_col),  
                subplots_positions = plots_positions,
                subplots_titles = filter_suptitle,
                subplots_save_ind = subplots_save_ind,
                legend_plot=False,
                x_data = x_sub_data,
                y_data = y_sub_data,
                axis_font_size = 26,
                fig_name = fig_name,
                fig_title = fig_title,
                fig_save_path = plot_path,
                x_axis_label = [('Redshift','Emitted Wavelenght [Å]')],
                x_ticks_dec = 2,
                x_axis_sides = 'both',
                x_top_axis_dec = 1,
                x_diff_fuct = ('RtW',filter_wl),
                y_axis_label = [(left_mag_prop['label'],right_mag_prop['label'])],
                y_axis_diff = True,
                y_axis_sides = 'both',
                y_axis_color = True,
                line_color=[(left_mag_prop['color'],right_mag_prop['color'])],
                line_alpha=[(left_mag_prop['alpha'],right_mag_prop['alpha'])],
                mark_style = [(left_mag_prop['marker'],right_mag_prop['marker'])],
                mark_size = [(left_mag_prop['marker_size'],right_mag_prop['marker_size'])],
                style_plot=[(left_mag_prop['plot_style'],right_mag_prop['plot_style'])],
                label_list = [(left_mag_prop['leg_lab'],right_mag_prop['leg_lab'])],
                zorder=[(left_mag_prop['zorder'],right_mag_prop['zorder'])],
                plot_show=False)



def plots_mag(system_list_mode,
              filters_list_mode,
              plot_max_columns,
              datacube_name_general,
              analysis_programme,
              galaxy_name,
              plot_save_mode,
              plot_path,
              filter_dict = {},
              filters_to_plot = [],
              mag_pairs_plot = [],
              wl_dict = {},
              plot_prop_dict = {}):     

    if system_list_mode == 0:
        filter_system_plot = filter_dict.keys()
    else:
        filter_system_plot = system_list_mode

    del_dict = copy.deepcopy(filter_dict)
    # For every filter system there is a plot sheet
    for filter_system in filter_system_plot:

        if filters_list_mode != 0:
            # Crear una lista con las claves a eliminar
            keys_a_eliminar = [key for key in del_dict[filter_system] if key not in filters_list_mode]
            # Eliminar las claves fuera del bucle
            for key in keys_a_eliminar:
                del del_dict[filter_system][key]

        # Obtener las claves actualizadas después de la eliminación
        filters_to_plot = list(del_dict[filter_system].keys())
        if len(del_dict[filter_system].keys()) == 0:
            filters_to_plot = keys_a_eliminar
        
        # For storing the plot positions
        plots_positions = []
        
        # For storing the plots in the correct
        # order to plot
        plot_order = []

        # For storing the filters names as 
        # titles of plots
        filter_suptitle = []
        
        # For storing the filters wavelenght
        # for convert from redshift to emitted wl
        filter_wl = []
        
        # Total columns of the plot, max is 6
        total_col  = 0
        
        # For store data in the correct order
        x_sub_data = []
        y_sub_data = []
        
        # Initializing the rows values
        # One row per filters system except if cols>6
        # for a given system
        rows = 0
        
        
        x_labels = []
        mag_plot_style = []
        mag_leg_label = []
        mag_labels = []
        mag_colors = []
        mag_markers = []
        mag_marker_sizes = []
        mag_alphas = []
        mag_zorder = []
        
        # For each pair of magnitudes in the same
        # sheet of plots
        for mag_pair in mag_pairs_plot:
            
            # Obtaining plot properties for each magntiude in
            # the correct format
            x_labels.append(('Redshift','Emitted Wavelenght [Å]'))
            mag_labels.append((plot_prop_dict[mag_pair[0]]['label'],
                               plot_prop_dict[mag_pair[1]]['label']))
            mag_plot_style.append((plot_prop_dict[mag_pair[0]]['plot_style'],
                               plot_prop_dict[mag_pair[1]]['plot_style']))
            mag_leg_label.append((plot_prop_dict[mag_pair[0]]['leg_lab'],
                               plot_prop_dict[mag_pair[1]]['leg_lab']))
            mag_colors.append((plot_prop_dict[mag_pair[0]]['color'],
                               plot_prop_dict[mag_pair[1]]['color']))
            mag_markers.append((plot_prop_dict[mag_pair[0]]['marker'],
                               plot_prop_dict[mag_pair[1]]['marker']))
            mag_marker_sizes.append((plot_prop_dict[mag_pair[0]]['marker_size'],
                               plot_prop_dict[mag_pair[1]]['marker_size']))
            mag_alphas.append((plot_prop_dict[mag_pair[0]]['alpha'],
                               plot_prop_dict[mag_pair[1]]['alpha']))
            mag_zorder.append((plot_prop_dict[mag_pair[0]]['zorder'],
                               plot_prop_dict[mag_pair[1]]['zorder']))
            
            
            # To avoiding duplicate filters
            filter_set = set()
            
            # Initializing column values
            cols = 0
            total_col_max = 0

                
            for filt in filters_to_plot:
                
                # If filter is in the list of filters to plot
                if filt in filters_to_plot and filt not in filter_set:
                    
                    # Avoiding to plot the same filter twice
                    filter_set.add(filt)
                        
                    # Name of each filter
                    filter_suptitle.append(filt)
                    
                    # Filter wavelength in the same order as the filters
                    filter_wl.append(wl_dict[filter_system][filt])
            
                    # Determining the order and the positions of the plots
                    plot_order.append(((filter_system,filt),(rows,cols)))
                    
                    plots_positions.append((rows,cols))
                    
                    # Obtainig the data for the plot
                    x_data = filter_dict[filter_system][filt]['z']
                    x_sub_data.append(x_data)
                    
                    # For the left Y axis
                    y_data_left = filter_dict[filter_system][filt][mag_pair[0]]
                    
                    # For the right Y axis
                    y_data_right = filter_dict[filter_system][filt][mag_pair[1]]
                    
                    # Tuple of y data values added to the y values list
                    y_sub_data.append((y_data_left,y_data_right))
                    
                    
                    # Determining the number of columns for the plot
                    total_col_max += 1
                    if total_col_max > total_col:
                        total_col = total_col_max
                    
                    # Correction for systems with more than 
                    # max_col filters to plot by adding one extra row
                    if cols // plot_max_columns != 0 and cols == plot_max_columns:
                        rows += 1
                        cols = 0
                        plots_positions.pop()
                        plots_positions.append((rows,cols))
                    else:
                        cols += 1
            
            # New row for each filter system        
            rows += 1
        
        # Total rows for the plot sheet
        total_row = rows
        
        fig_name = f'{datacube_name_general}_{filter_system}'
        fig_title = (f'\\textbf{{{analysis_programme}}} output for galaxy \\textbf{{{galaxy_name}}} '
                     f'with \\textbf{{{filter_system}}} filter system')
        
        if len(filter_set) == 0:
            break
        
        elif len(filter_set) == 1 or plot_save_mode==1:
            fig_title = (f'\\textbf{{{analysis_programme}}} output for galaxy \\textbf{{{galaxy_name}}} \n'
                     f'with \\textbf{{{filter_system}}} filter system')
        
        
        # Plotting the filters with this configuration
        plot_gen(subplots_conf = True,
                subplots_num=(total_row,total_col),
                subplots_positions = plots_positions,
                subplots_titles = filter_suptitle,
                #subplots_split_groups = len(filter_dict.keys()),
                subplots_split_groups = len(mag_labels),
                legend_plot=False,
                x_data = x_sub_data,
                y_data = y_sub_data,
                fig_fonsize = 20,
                axis_font_size = 26,
                fig_title = fig_title,
                fig_name = fig_name,
                fig_save_path = plot_path,
                aspect_ratio=[1,1],
                x_axis_label = x_labels,
                x_ticks_dec = 2,
                x_axis_sides = 'both',
                x_top_axis_dec = 1,
                x_diff_fuct = ('RtW',filter_wl),
                y_axis_label = mag_labels,
                y_axis_diff = True,
                y_axis_sides = 'both',
                y_axis_color = True,
                line_color=mag_colors,
                line_alpha=mag_alphas,
                mark_style = mag_markers,
                mark_size = mag_marker_sizes,
                style_plot=mag_plot_style,
                label_list=mag_leg_label,
                zorder=mag_zorder,
                plot_show=False)           
        
        

def run_plotting_pipeline():
    
    cronometro = Cronometro()
    #cronometro.iniciar()
    
    # Obtaingin configuration values
    (analysis_programme,system_list_mode,filters_list_mode,
     mag_pairs_mode,plot_mode,plot_save_mode,plot_max_columns,
     plot_analysis_version) = argparse_values(phase='plot')

    # Obtaining the current working directory
    # in which the programme is executed
    SCRIPT_PATH = Path(__file__).resolve()
    cwd_scripts = SCRIPT_PATH.parent
    cwd_src = cwd_scripts.parent
    cwd_global = cwd_src.parent
    
    datacube_folder_name_list = []
    datacube_folder_path_list = []
    
    for path, subdir, file in sorted(os.walk(cwd_global)):
        for subdir_name in subdir:
            if 'original_fits' in subdir_name:
                datacube_folder_name_list.append(subdir_name)
                datacube_folder_path_list.append(path)
    
    for path_pos,path in enumerate(datacube_folder_path_list): 
        
        # Changing the directory to obtain the original path
        os.chdir(path)
        cwd_galaxy = os.getcwd()
        
        # Returning to the scripts directory
        os.chdir(cwd_scripts)

        # Galaxy name obtaiane from datacube directories
        galaxy_name = (datacube_folder_name_list[path_pos].split('/')[-1]).split('_')[0]
        
        # Absolute path to the original fits datacubes
        # This files are the ones to be analyzed
        datacube_original_folder_path = f'{cwd_galaxy}/{galaxy_name}_original_fits'

        datacube_name_list = []
        for datacube_name in sorted(os.listdir(datacube_original_folder_path)):
            if '.fits' in datacube_name:
                datacube_name_list.append(datacube_name)

        # Datacubes general name for the galaxy without the filter name
        datacube_name_general = '_'.join((datacube_name_list[0].split('.fits')[0]).split('_')[:-1])
        
        version_file_name = f'txt_{galaxy_name}_version.txt'
        #version_file_path = f'{cwd_global}/{galaxy_name}/{version_file_name}'
        version_file_path = f'{cwd_global}/data/{version_file_name}'
        
        # Obtaining the latest version runned
        if plot_analysis_version == 'latest':
            last_version = version_file_last(version_file_path)

        else:
            last_version = plot_analysis_version
            
        if last_version < 10:
                last_version = f'0{last_version}'
                
        analysis_folder_name = f'{galaxy_name}_analysis_V{last_version}_{analysis_programme}'
        analysis_folder_path = f'{cwd_galaxy}/{analysis_folder_name}'
            
        # Path to the csv with the data
        csv_folder_name = f'{analysis_folder_name}_csv'
        csv_folder_path = f'{analysis_folder_path}/{csv_folder_name}'
        
        datacube_csv_filter_list = []
        for datacube_name in sorted(os.listdir(csv_folder_path)):
            if '.csv' in datacube_name:
                filt_or_name = datacube_name.split('.fits')[0].split('_')[-2]
                datacube_csv_filter_list.append(filt_or_name)

        # Path to the folder of the plots
        plot_folder_name = f'{analysis_folder_name}_plots'
        plot_path = f'{analysis_folder_path}/{plot_folder_name}'
        create_folder(plot_path,overwrite=False)
            
        # The names of the filter systems in alphabetic order
        # Dictionary with the corresponding wavelength of each filter
        # http://svo2.cab.inta-csic.es/theory/fps/               
        filter_system_names,filter_wavelength_dict,filter_names_dict = filter_wl_dict()    
        
        
        # Obtaining the filters sorted in a dictionary and grouped by 
        # filter system
        (df_filter_dict_short,df_filter_dict_short_wl,
        df_filter_dict_large_to_short) = filter_sort_dict(csv_folder_path,
                                                          filter_wavelength_dict,
                                                          filter_names_dict)
        
        
        # Plot selected filters
        if filters_list_mode != 0 and isinstance(filters_list_mode, list)==True:
        
            # Filters to plot
            filters_to_plot = filters_list_mode
        
        # Plot all available filters
        elif filters_list_mode == 0:
            filters_to_plot = []
            for filt_name in datacube_csv_filter_list:
                filters_to_plot.append(df_filter_dict_large_to_short[filt_name])
        
        
        # Plot selected magnitudes
        if mag_pairs_mode != 0 and isinstance(mag_pairs_mode, list)==True:
        
            # Magnitudes to plot
            mag_pairs_plot = [('Ind','Eff_rad_kpc'),('Mag','Mag'),('Ax_rat','Pos_ang')]
        
        # Plot all available filters
        elif mag_pairs_mode == 0:
            if analysis_programme == 'Galfit':
                # Magnitudes to plot
                mag_pairs_plot = [('Ind','Eff_rad_arcsec'),('Mag','Mag'),('Ax_rat','Pos_ang'),('X_cent','Y_cent')]
            elif analysis_programme == 'Imfit':
                mag_pairs_plot = [('Ind','Eff_rad_arcsec'),('I_e','I_e'),('Ax_rat','Pos_ang'),('X_cent','Y_cent')]
        
        # Obtaining plot properties for each magnitude
        mag_plot_prop_dict = plot_prop_conf(mag_pairs_plot)
        
        # Plot mode
        try:
            if plot_mode == 0 or plot_mode == 2 or plot_save_mode == 1:
                # If all filters have to be plotted in the same sheet
                plots_filter(system_list_mode,
                             filters_list_mode,
                             datacube_name_general,
                             analysis_programme,
                             galaxy_name,
                             plot_save_mode,
                             plot_max_columns,
                             plot_path,
                             filter_dict=df_filter_dict_short, 
                             filters_to_plot=filters_to_plot,
                             mag_pairs_plot = mag_pairs_plot,
                             wl_dict=df_filter_dict_short_wl,
                             plot_prop_dict=mag_plot_prop_dict)
            
            if (plot_mode == 0 or plot_mode == 2) and plot_save_mode == 2:
                plot_save_mode = 0
                # If all filters have to be plotted in the same sheet
                plots_filter(system_list_mode,
                             filters_list_mode,
                             datacube_name_general,
                             analysis_programme,
                             galaxy_name,
                             plot_save_mode,
                             plot_max_columns,
                             plot_path,
                             filter_dict=df_filter_dict_short, 
                             filters_to_plot=filters_to_plot,
                             mag_pairs_plot = mag_pairs_plot,
                             wl_dict=df_filter_dict_short_wl,
                             plot_prop_dict=mag_plot_prop_dict)
            
            if plot_mode == 1 or plot_mode == 2:
                # If there is a sheet for each Telescope filter System
                plots_mag(system_list_mode,
                          filters_list_mode,
                          plot_max_columns,
                          datacube_name_general,
                          analysis_programme,
                          galaxy_name,
                          plot_save_mode,
                          plot_path,
                          filter_dict=df_filter_dict_short,
                          filters_to_plot=filters_to_plot,
                          mag_pairs_plot = mag_pairs_plot,
                          wl_dict=df_filter_dict_short_wl,
                          plot_prop_dict=mag_plot_prop_dict)
            
            print('Plots are done')

        except Exception as e:
            print(f'\nError plotting results:\n{e}')