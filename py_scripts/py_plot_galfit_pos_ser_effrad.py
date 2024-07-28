
# This script creates the plots for the positon of galaxy center
# for the Sersic index and the effective radius

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FixedLocator
import functools
from itertools import groupby

import os
import subprocess
import shutil
import re

import pdb
import time

from py_galfit_filter_wavelenght import filter_wl_dict


#-------------------FUNCTIONS-------------------#

# FUNCTIONS FOR THE SECONDARY AXES TICKS
# Transforming the redshift into wavelenght
def RtW(z,filter_wl):

    wavelenght = (filter_wl/10)/(1+z)
    
    np.seterr(divide = 'ignore') 
    np.seterr(invalid = 'ignore')
    
    return wavelenght
    
# Transforming the wavelenght into redshift
def WtR(wl,filter_wl):
    
    redshift = ((filter_wl/10)/wl) - 1 

    np.seterr(divide = 'ignore') 
    np.seterr(invalid = 'ignore')
    
    return redshift

# Position of the minor ticks in top axis
def minor_ticks_pos(z_labels,filter_value):
    tick_pos_list = []
    
    for i in range(len(z_labels)-1):
        tick_pos_z = np.linspace(z_labels[i],z_labels[i+1],5)
        for i in range(len(tick_pos_z)-1):
            tick_pos_list.append(RtW(tick_pos_z[i],filter_value))

    return tick_pos_list

# PLOT LAYOUT FUNCTIONS
# Creating the grid layour for the array plot    
def layout(row,col):

    if col<5 or col==5:
        
        return row,col
        
    else:
        
        return row*2,(col//2+col%2)
    

# Sorting the plots posititons in the correct way
def grid_pos(row,col):

    positions = []
    
    # obtaining the grid positions
    for i in range(col):
        for j in range(row):
            positions.append((j,i))
    
    # if there are only 2 rows do nothing
    if row<2 or row == 2:
    
        return positions
        
    # if there are more that 2 rows we want at first to 
    # put the plots in the row 1 and 3 and then move
    # to another column. Then put the plots in the 2 and 4 rows
    else:
        even_pos = []
        odd_pos = []
        
        # the 1 and 3 rows correspond to the even positions tuples
        # from the grid positions list
        for i  in range(len(positions)):
            if i%2 == 0:
                even_pos.append(positions[i])
            
            else: 
                odd_pos.append(positions[i])
                
        pos_sort = even_pos+odd_pos
        
        return pos_sort


# OTHER FUNCTIONS
# Obtaining the key of a dictionary from a value
def get_key(val,dic):
    
    for key, value in dic.items():

        if val == value:
            return key



#-------------------MAIN CODE-------------------#

# Obtaining the current working directory and the galaxy
cwd = os.getcwd()
galaxy = cwd.split('/')[-1].split('_')[0]
psf = cwd.split('/')[-1].split('_')[2]
version = cwd.split('/')[-1].split('_')[-1]

# The names of the filter systems in alphabetic order
# Dictionary with the corresponding wavelenght of each filter
# http://svo2.cab.inta-csic.es/theory/fps/               
filter_system_names,filter_wavelenght_dict,filter_names_dict = filter_wl_dict()        

# lists to store the filters names and the csv paths 
filter_name_list = []
file_path_list = []

# walking thorught the tree directory
for (dirpath, dirnames, filenames) in os.walk(cwd):

    # sorting the filters directories in alhpabetical order
    dirnames.sort()
    for subdir in dirnames:
        
        for file in filenames:
            
            if '.csv' in file:
                
                fits_name = dirpath.split('/')[-1]
                
                filter_name = fits_name.split('_')[-1]
                
                if filter_name[-1] == 'M':
                    filter_name = filter_name.replace('M','W')
                    
                if filter_name[:3] == 'HST':
                    filter_name = filter_name.replace('W','H')
                
                filter_name_list.append(filter_name)
                            
                file_path = os.path.join(dirpath,file)
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
                'HSTF475W', 'HSTF555W', 'HSTF569W', 'HSTF606W', 'HSTF791W', 'HSTF814W')):
            filter_name_group[i][j] = f'{flt[:-1]}W'            
            
# Sorting the system group wavlenghts from min to max
filter_wavelenght_group = [[] for i in range(len(filter_name_group))]
filter_wavelenght_group_title = [[] for i in range(len(filter_name_group))]

for index_system,filter_system in enumerate(filter_name_group):
    for index_filtr,filtr in enumerate(filter_system):
        
        filter_wavelenght_group[index_system].append(filter_wavelenght_dict[filtr])
        filter_wavelenght_group_title[index_system].append(filter_names_dict[filtr])

        filter_wavelenght_group[index_system] = sorted(
                            filter_wavelenght_group[index_system])
        filter_wavelenght_group_title[index_system] = sorted(
                            filter_wavelenght_group_title[index_system])

# Sorting the filters name in increasing wavelenght
filter_name_group_sorted = [[] for i in range(len(filter_name_group))]
filter_name_group_title_sorted = [[] for i in range(len(filter_name_group))]

for index_wl_system,wl_system in enumerate(filter_wavelenght_group):
    
    for index_wl_filtr,wl_filter in enumerate(wl_system):
        
        key = get_key(wl_filter,filter_wavelenght_dict)
        key_title = filter_names_dict[key]

        filter_name_group_sorted[index_wl_system].append(key)
        filter_name_group_title_sorted[index_wl_system].append(key_title)

# Creating an array plot for every filter system
for filter_system_index, filter_system in enumerate(filter_name_group_sorted):
    
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
    df_filter = []    

    for filter_name in filter_system:

        # Obtaining the csv path
        path_csv = [path for path in file_path_list if filter_name in path][0]
        
        # Storing the dataframes for each csv in a list
        df_filter.append((str(filter_system[df_index]),pd.read_csv(path_csv, sep=',')))
        
        df_index = df_index + 1
    
    #-------------------PLOTTING THE RESULTS-------------------#
    
    # We are going to create a plot array with the center positions,
    # the sersic index and the effective radius for each filter system
    
    print(f'Creating an array of plots for {galaxy} with {filter_name_title} filter system\n')
    
    # Using tuples because we have a pair of x-axes and a apair of y-axes
    
    # The second z will be used to compute the corresponding wavelenght for the 
    # secondary x-axis labels
    x_axis_magnitude = [('z','z'),('z','z')]
    x_axis_label = [('Redshift','Wavelenght Emitted [nm]'),
                    ('Redshift','Wavelenght Emitted [nm]')]
    
    # We can define the pair of magnitudes we want to plot togheter on y-axis
    y_axis_magnitude = [('X_cent','Y_cent'),('Ind','Eff_rad_kpc')] 
    y_axis_label = [('X center position [pix]','Y center position [pix]'),
                    ('Sersic Index','Effective radius [Kpc]')] 
    
    
    # calling the layout func to define the correct grid shape
    total_row,total_col = layout(len(x_axis_magnitude),len(filter_system))
    
    # Defining an array layout for the plots (rows, columns)
    # with a corresponding size
    fig, axs = plt.subplots(nrows=total_row, ncols=total_col,
                figsize=(total_col*5, total_row*5))
    
    # Calling the grid_pos func to define the correct position order
    positions = grid_pos(total_row,total_col)
    
    # Defining and index for the positions
    pos_index = 0
    
    # Title of the plot
    fig.suptitle(f'{galaxy} with {filter_name_title} filters system with the {psf} psf {version}',
            ha='center', va='top', fontweight='bold', fontsize='16')
    
    if len(filter_system) != len(positions)/2:
        
        for i in positions[-2:]:
            fig.delaxes(axs[i[0]][i[1]])
        
    # Number of filters on each filter system
    if len(filter_system) > 0:
    
        for filters in range(len(filter_system)):
            
            # Number of magnitudes we want to plot
            for mag in range(len(x_axis_magnitude)):

                # Name of the current plotted filter
                plot_filter = filter_name_group_title_sorted[filter_system_index][filters]

                # Defining short variables for x magnitudes from bottom x-axis
                x_mag_1 = x_axis_magnitude[mag][0]
                x_mag_2 = x_axis_magnitude[mag][1]
                
                x_lab_1 = x_axis_label[mag][0]
                x_lab_2 = x_axis_label[mag][1]
                
                # Defining short variables for y magnitudes from left y-axis
                y_mag_1 = y_axis_magnitude[mag][0]
                y_mag_2 = y_axis_magnitude[mag][1]
                
                y_lab_1 = y_axis_label[mag][0]
                y_lab_2 = y_axis_label[mag][1]
                
                # Bottom x-axis values
                x_values = (df_filter[filters][1])[f'{x_mag_1}'].to_numpy()

                # Create the subplot
                if len(axs.shape)>1:
                    ax = axs[positions[pos_index][0],positions[pos_index][1]]
                else:
                    ax = axs[positions[pos_index][0]]

                # Adding the filter name to every subplot
                ax.set_title(plot_filter, fontweight='bold', fontsize='13')
                
                # Left y-axis values
                y_values = (df_filter[filters][1])[f'{y_mag_1}'].to_numpy()
                y_values_delta = (df_filter[filters][1])[f'{y_mag_1}_err'].to_numpy()        
                   
                # Color of data
                if y_mag_1 == 'X_cent':
                    color = 'darkorange'

                else:
                    color = 'red'

                # Plot data
                ax.scatter(x_values, y_values, color=color, marker='.',s=20)
                
                # Customize bottom x-axis and left y-axis with labels
                ax.set_xlabel(f'{x_lab_1}', fontweight='normal', fontsize='13')
                ax.set_ylabel(f'{y_lab_1}', color=color, fontweight='bold', fontsize='13')
                ax.tick_params(axis='y', labelcolor=color)
                ax.tick_params(axis='both', which='major', labelsize=11)
                
                # Number of ticks in the axex and space between
                # the extreme values and the edges of the plot            
                ax.set_xticks(np.linspace(np.min(x_values), np.max(x_values), 5))
                ax.set_xmargin(0.1)

                
                if all(item == y_values[0] for item in y_values) == False:

                    ax.set_yticks(np.linspace(np.min(y_values), np.max(y_values), 5))
                    
                else: 
                
                    ax.set_yticks(np.linspace(np.min(y_values), np.max(y_values)+1,1))
                    
                
                ax.set_ymargin(0.05)
                
                # Allowing minor ticks without text
                ax.minorticks_on()
                
                # Number of decimals in the ticks
                ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                
                # Grid in horizontal and vertical orientations
                ax.grid(True)
                
                # Top x-axis and its ticks
                # https://stackoverflow.com/questions/68715304/dual-x-axis-with-same-data-different-scale
                # https://docs.python.org/3.8/library/functools.html#functools.partial
                
                # Obtainin the wavelenght of the corresponding filter
                filter_wl_value = filter_wavelenght_dict[filter_system[filters]]
                
                # Creating the secondary top x-axis with
                # the ticks in wavelenght using the partial method to pass
                # arguments to the RtW and WtR functions
                axxtop = ax.secondary_xaxis('top', 
                    functions=(functools.partial(RtW, filter_wl=filter_wl_value),
                               functools.partial(WtR, filter_wl=filter_wl_value)))
                               
                # From the bottom x-axis ticks in redshift we can compute the 
                # same values but in wavelenght with RtW function
                z_ticks = ax.get_xticks()
                wl_ticks = RtW(z_ticks,filter_wl_value)
                
                axxtop.set_xticks(wl_ticks)

                # Minor ticks positions of top axis
                ticks_pos = minor_ticks_pos(z_ticks,filter_wl_value)
                axxtop.xaxis.set_minor_locator(FixedLocator(ticks_pos))
                
                # Allowing minor ticks without text
                axxtop.minorticks_on()
                       
                # Label for this axis        
                axxtop.set_xlabel(f'{x_lab_2}', fontweight='normal', fontsize='13')
                axxtop.tick_params(axis='x', which='major', labelsize=11)
                
                # Number of decimals in the ticks
                axxtop.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                
                
                # Right y-axis values for the second magnitude to plot
                y_values = (df_filter[filters][1])[f'{y_mag_2}'].to_numpy()
                y_values_delta = (df_filter[filters][1])[f'{y_mag_2}_err'].to_numpy()    
                
                
                # Creating a secondary y-axis in the right
                # that shares the x-axis with the left one
                axyrig = ax.twinx() 
                
                # Color of data
                if y_mag_2 == 'Y_cent':
                    color = 'limegreen'
                else:
                    color = 'blueviolet'

                # Plot data
                axyrig.scatter(x_values, y_values, color=color, marker='x',s=10)
                
                # Customize right y-axis with label
                axyrig.set_ylabel(f'{y_lab_2}', color=color, fontweight='bold', fontsize='12')
                axyrig.tick_params(axis='y', labelcolor=color)
                axyrig.tick_params(axis='y', which='major', labelsize=11)
                
                # Number of ticks in the axis, the space between the 
                # data and the edges and the decimals in the ticks
                if all(item == y_values[0] for item in y_values) == False:

                    axyrig.set_yticks(np.linspace(np.min(y_values), np.max(y_values), 5))
                    
                else: 
                
                    axyrig.set_yticks(np.linspace(np.min(y_values), np.max(y_values)+1,1))
                
                axyrig.set_ymargin(0.05)
                axyrig.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                
                axyrig.grid()
                
                # Allowing minor ticks without text
                axyrig.minorticks_on()
                    
                # Moving to the next position in the grid layout
                pos_index = pos_index + 1
             
            
    # Adding enought space for the plots and the title
    fig.tight_layout(pad=2)

    # Saving the plots in pdf
    plt.savefig(f'{cwd}/{galaxy}_plot_galfit_{filter_name_title}_{psf}_{version}_pos_ser_radeff.pdf', format='pdf', dpi=1000, bbox_inches='tight')

    # Closing the plot to avoid overlapse
    plt.close()    


print('All galfit plots are done')
