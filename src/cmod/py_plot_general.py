### Imports:
#############################################################################################
import sys

import pdb

import numpy as np

import pandas as pd
pd.set_option('mode.chained_assignment', None)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FixedLocator
import matplotlib.font_manager as font_manager
import functools


from py_convert_functions import RtW, WtR


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