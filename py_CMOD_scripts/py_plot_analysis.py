import pandas as pd
from itertools import groupby
import sys
import argparse
import os
import pdb
import copy


from py_galfit_filter_wavelength import filter_wl_dict
from py_create_folder import create_folder
from py_plot_general import plot_gen
from py_plot_configration import plot_prop_conf
from py_cronometro import Cronometro
from py_config_argparse import argparse_values
from py_version_control import version_directory, version_file, version_file_last

#-------------------FUNCTIONS-------------------#

# Obtaining the key of a dictionary from a value
def get_key(val,dic):
    
    for key, value in dic.items():
        if val == value:
            return key

#-------------------MAIN CODE-------------------#

def filter_sort_dict():
    
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
  
def plots_filter(filter_dict = {},
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


def plots_mag(filter_dict = {},
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
            
        pdb.set_trace()
        
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





if __name__ == '__main__':
    
    cronometro = Cronometro()
    #cronometro.iniciar()
    
    # Obtaingin configuration values
    (analysis_programme,system_list_mode,filters_list_mode,
     mag_pairs_mode,plot_mode,plot_save_mode,plot_max_columns,
     plot_analysis_version) = argparse_values(phase='plot')

    # Obtaining the current working directory
    # in which the programme is executed
    cwd_scripts = os.getcwd()
    cwd_global = cwd_scripts + '/..'
    
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
        version_file_path = f'{cwd_global}/{galaxy_name}/{version_file_name}'
        
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
        df_filter_dict_large_to_short) = filter_sort_dict()
        
        
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
                plots_filter(filter_dict=df_filter_dict_short, 
                            filters_to_plot=filters_to_plot,
                            mag_pairs_plot = mag_pairs_plot,
                            wl_dict=df_filter_dict_short_wl,
                            plot_prop_dict=mag_plot_prop_dict)
            
            if (plot_mode == 0 or plot_mode == 2) and plot_save_mode == 2:
                plot_save_mode = 0
                # If all filters have to be plotted in the same sheet
                plots_filter(filter_dict=df_filter_dict_short, 
                            filters_to_plot=filters_to_plot,
                            mag_pairs_plot = mag_pairs_plot,
                            wl_dict=df_filter_dict_short_wl,
                            plot_prop_dict=mag_plot_prop_dict)
            
            if plot_mode == 1 or plot_mode == 2:
                # If there is a sheet for each Telescope filter System
                plots_mag(filter_dict=df_filter_dict_short,
                        filters_to_plot=filters_to_plot,
                        mag_pairs_plot = mag_pairs_plot,
                        wl_dict=df_filter_dict_short_wl,
                        plot_prop_dict=mag_plot_prop_dict)
            
            print('Plots are done')
        
        except Exception as e:
            print(f'\nError plotting results:\n{e}')