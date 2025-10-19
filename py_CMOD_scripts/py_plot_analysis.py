import pandas as pd
from itertools import groupby
import os
from pathlib import Path

from cmod.io import argparse_values
from cmod.photometry import filter_wl_dict
from cmod.plotting import plot_prop_conf, filter_sort_dict, plots_filter, plots_mag
from cmod.utils import Cronometro, create_folder, version_file_last

  

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
    cwd_scripts = SCRIPT_PATH.parent # py_CMOD_scripts/
    cwd_global = cwd_scripts.parent  # CMOD/
    DATA_DIR = cwd_global / "data"   # Define data dir
    
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

if __name__ == '__main__':
    try:
        run_plotting_pipeline()
    except Exception as e:
        print(f'\nError running plotting pipeline:\n{e}')
    
        
