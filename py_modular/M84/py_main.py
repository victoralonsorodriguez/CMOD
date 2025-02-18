
import numpy as np
import pandas as pd
import scipy as scp

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import functools

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import os
import re
import sys
import subprocess
import shutil
import mmap
import struct

import pdb
import time
import warnings
warnings.filterwarnings('ignore')

# Packages for initial condition
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.signal import deconvolve

from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture
from photutils.isophote import Ellipse

from py_galfit_constraints_script import create_constraints
from py_galfit_constraints_values import constraints_values
from py_galfit_initial_params import initial_params
from py_galfit_initial_params_script import create_initiaL_params
from py_galfit_script import create_script
from py_galfit_kpc_correction import kpc_correction

from py_plot_general import plot_gen
from py_open_fits import open_fits
from py_create_folder import create_folder
from py_mag_counts_convert import fits_mag_to_counts, fits_counts_to_mag, values_counts_to_mag
from py_rad_deg import deg_to_rad,rad_to_deg,rad_to_deg_abs
from py_round_number import round_number
from py_max_center_value import max_center_value
from py_dataframe import galfit_init_dataframe, galfit_create_dataframe
from py_isophote_fitting import isophote_fitting
from py_initial_condition import initial_conditions


def main():
    
    '''
    This function contains the analysis process
    '''
    
    # The analysis is carried out for each datacube
    for datacube_pos,datacube_name in enumerate(datacube_name_list):
        
        datacube_name_noext = datacube_name.split('.fits')[0]
        
        # Creating a directory with the datacube name to store the 
        # files produced during the analysis
        datacube_folder_path = f'{cwd}/{datacube_name_noext}'
        create_folder(datacube_folder_path,
                      overwrite=True)
        
        # Copying the datacube from the original folder
        # to its corresponding folder
        datacube_original_path = f'{datacube_original_folder_path}/{datacube_name}'
        shutil.copyfile(datacube_original_path,
                        f'{datacube_folder_path}/{datacube_name}')

        # The datacube is in magnitudes but for the analysis is
        # required its value in counts       
        # The images have a padding that should be removed 
        datacube_counts_path =  fits_mag_to_counts(fits_input_path = datacube_original_path,
                                                   fits_output_path = datacube_folder_path,
                                                   const = pxsc_zcal_const,
                                                   remove_padding = True)
        
        # We need to open the fits to analyse each frame
        datacube_flux_hdr, datacube_flux_data, datacube_flux_name  = open_fits(datacube_counts_path)
        
        # Redshift step of the datacube
        redshift_step = datacube_flux_hdr['CDELT3']
        
        
        #--------CREATING A DATAFRAME TO STORE THE DATA------#

        df = galfit_init_dataframe()

        # Export dataframe to text file
        csv_file = open(f'{cwd}/{datacube_name_noext}_galfit.csv','w+')
        
        # We should analyse each frame separatelly
        
        analysis_range = datacube_flux_data.shape[0]
        analysis_range = 3
        
        for frame_pos in range(analysis_range):
            
            # Obtaining the corresponding redshift of each frame
            redshift_value = round_number(frame_pos * redshift_step,2)
            
            # We are going to export each frame of the data cube individually
            frame = datacube_flux_data[frame_pos]
            
            frame_name = f'{datacube_name.split(".")[0]}_z{redshift_value:.2f}_counts.fits'
            frame_name_noext = frame_name.split('.fits')[0]
            frame_path = f'{cwd}/{frame_name}'
                        
            # Depending on the used version we are going to use a median 
            # filter in order to reduce the image noise
            # By default the filter is aplied
            if ('--medfilt_inout' in sys.argv or '--medfilt' in sys.argv or
                '--inout' not in sys.argv or '--original' not in sys.argv):
                
                # Median filter applied to each frame
                frame = scp.signal.medfilt(frame, kernel_size = 3)

            # saving the subframe
            fits.writeto(f'{frame_path}', frame, header=datacube_flux_hdr, overwrite=True)
            
            
            # Obtaining the center of the galaxy as the maximun central pixel value
            # Only for the first image
            max_pix_value_center, max_pix_value_center_pos = max_center_value(frame,crop_factor)
            
            '''
            Initial conditions for each datacube
            '''
            
            if frame_pos == 0:
                
                #gal_iso_fit_csv_path = isophote_fitting(frame_path,max_pix_value_center_pos,cons=pxsc_zcal_const)
                gal_iso_fit_csv_path = './m84_VBIN018_SL_zSimJ_EucHab_z0.00_counts_isophote.csv'
                gal_iso_fit_df = pd.read_csv(gal_iso_fit_csv_path)
                break_pos, break_pos_bul, I_0_disk_in, h_disk_in, n_bul_in, r_e_bul_in, I_e_bul_in = initial_conditions(frame_name_noext,
                                                                                                                        gal_iso_fit_df, 
                                                                                                                        'sma', 
                                                                                                                        'intens', 
                                                                                                                        'intens_err',
                                                                                                                        const=pxsc_zcal_const)
                ell_mean_in = gal_iso_fit_df.loc[:, 'ellipticity'].mean()
                ax_rat_in = round_number(1 - ell_mean_in,3)
                pa_rad_mean_in = gal_iso_fit_df.loc[:, 'pa'].mean()
                pa_deg_in = round_number(rad_to_deg_abs(pa_rad_mean_in)-90,3)
            
            
            # Obtaining the frame shape
            frame_shape = frame.shape
            
            # Obtaining them toal surface magnitude of the frame
            frame_total_flux = np.nansum(frame)
            frame_total_mag = values_counts_to_mag(frame_total_flux,
                                                    const = pxsc_zcal_const)    
            
            # Creating a script to run galfit
            script_name = f'{frame_name_noext}.script'
            script_path = f'{cwd}/{script_name}'
            script_file = open(f'{script_path}','w+')
            
            # the ouput file name will be
            output_file = f'{frame_name_noext}_galfit.fits'
            output_file_path = f'{cwd}/{output_file}'
        
            
            # For versions using inout methodology
            if ('--inout' in sys.argv or '--medfilt_inout' in sys.argv or
                '--original' not in sys.argv or '--medfilt' not in sys.argv):
    
                # For all frames except for the first one we introduce the output 
                # parameters of previous frame as input
                if frame_pos != 0:

                    # calling the function to create the script
                    create_script(file=script_file,
                                input_file_name=frame_name,
                                output_file_name = output_file,
                                file_shape = frame_shape,
                                zp_const = zcal_const,
                                pix_scl = pixel_scale,
                                img_center = (y_center,x_center),
                                int_mag = mag,
                                eff_rad = eff_rad,
                                ser_ind = ser_index,
                                ax_rat = ax_rat,
                                pos_ang = pos_ang,
                                psf= psf_fits_name,
                                cons_file = constraints_file_name)

                # Initial parameters for the first frame of the filter datacube
                else:

                    # calling the function to create the script
                    create_script(file=script_file,
                                input_file_name=frame_name,
                                output_file_name = output_file,
                                file_shape = frame_shape,
                                zp_const = zcal_const,
                                pix_scl = pixel_scale,
                                img_center = max_pix_value_center_pos,
                                int_mag = frame_total_mag,
                                eff_rad = r_e_bul_in,
                                ser_ind = n_bul_in,
                                ax_rat = ax_rat_in,
                                pos_ang = pa_deg_in,
                                psf= psf_fits_name,
                                cons_file = constraints_file_name)

            # For versions not using inout methodology
            elif '--original' in sys.argv or '--medfilt' in sys.argv:

                # Initial parameters for every frame of the filter datacube
                create_script(file=script_file,
                            input_file_name=frame_name,
                            output_file_name = output_file,
                            file_shape = frame_shape,
                            zp_const = zcal_const,
                            pix_scl = pixel_scale,
                            img_center = max_pix_value_center_pos,
                            int_mag = frame_total_mag,
                            eff_rad = r_e_bul_in,
                            ser_ind = n_bul_in,
                            ax_rat = ax_rat_in,
                            pos_ang = pa_deg_in,
                            psf= psf_fits_name,
                            cons_file = constraints_file_name)

            
            # closing the script file
            script_file.close()

            # runing galfit for the subframe fits with the corresponding
            subprocess.run(['galfit', script_name])
            
        
            #-----OBTAINING THE PARAMETERS OF THE ANALYSIS----#
            
            # Extracting the galfit results data from the header of each output model frame
            
            # This 'if' is here because sometimes galfit doesn't fit a parameter and due that
            # it doesn't create an output file. This 'if' is to avoid the error of file not found
            if os.path.isfile(f'{output_file_path}'):
                
                df,y_center,x_center,mag,eff_rad,ser_index,ax_rat,pos_ang = galfit_create_dataframe(df,galaxy_name,output_file_path,pxsc_zcal_const)
                
        #----WRITING THE DATAFRAME INTO A CSV----#

        # Sort the dataframe by redshift
        df.sort_values(by=['z'],ascending=True, inplace=True, ignore_index=True)
        
        df_string = df.to_csv(header=True, index=False, sep=',')

        csv_file.write(df_string)

        csv_file.close()
        
        #----ORGANIZING ALL THE CREATED FILES INTO NEW FOLRDES----#
        
        # moving all the output and created files to a folder in order to keep the order in the directory        
        for output_file in sorted(os.listdir(cwd)):

            if ('.py' not in output_file and os.path.isdir(output_file) == False and
                ('galfit.' in output_file or 'fit.' in output_file or datacube_name_noext in output_file)):
                
                os.replace(f'{cwd}/{output_file}', f'{datacube_folder_path}/{output_file}')
                continue
                
            if '.csv' in output_file:

                shutil.copyfile(f'{cwd}/{output_file}',f'{datacube_folder_path}/{output_file}')


if __name__ == '__main__':
    
    '''
    This part of the code corresponds with global options

    The variables defined here can be used for all the functions
    '''
    
    # Computing the total analysis time
    start_time = time.time()
    
    # Obtaining the current working directory
    # in which the programme is executed
    cwd = os.getcwd()
    previus_cwd = cwd + '/..'
    
    # Obtaining the name of the galaxy
    galaxy_name = (cwd.split('/')[-1]).split('_')[0]
    
    # Absolute path to the original fits datacubes
    # This files are the ones to be analyzed
    datacube_original_folder_path = f'{cwd}/{galaxy_name}_original_fits'
    
    # Obtaining the datacube names
    datacube_name_list = []
    for datacube_name in sorted(os.listdir(datacube_original_folder_path)):
        if '.fits' in datacube_name:
            datacube_name_list.append(datacube_name)
            
    
    # We need the  pixel scale and the calibration constant to convert 
    # from counts to magnitudes
    pixel_scale = 0.2 # arcsec/pix
    zcal_const = 25
    pxsc_zcal_const = (pixel_scale,zcal_const)
    
    crop_factor = 10

    if '--version_control' in sys.argv:

        # Checking the latest version run
        version_file_path = f'{previus_cwd}/{galaxy_name}_version.txt'
        version_file = open(f'{version_file_path}','r')
        last_line_version_file = str(subprocess.check_output(['tail', '-1', version_file_path]))[2:-1]

        folder_new_version_name = last_line_version_file.split(' # ')[0]
        folder_new_version_path = f'{previus_cwd}/{folder_new_version_name}'

        # Depending on the folder name a large PSF or a
        # small psf will be created
        folder_name_elements = folder_new_version_name.split('_')
        version = folder_name_elements[1]

    # PSF: Point Spread Function 

    # Size of the PSF for galfit
    if 'psf' not in sys.argv:
    #if 'psf' not in folder_name_elements:
        psf_size = 'large'
    else:
        psf_size = 'small'
        
    folder_new_version_path = cwd
    version = 0


    # Name of the psf
    psf_fits_name = f'psf_{galaxy_name}_{psf_size}_flux_norm.fits'
    psf_fits_path = f'{folder_new_version_path}/{psf_fits_name}'

    # If an specific psf is not created for the galaxy
    # it will be used a general small
    if os.path.exists(f'{psf_fits_path}') == False:

        psf_fits_name = 'psf_small_stacked_original_flux_norm.fits'
        shutil.copyfile(f'./{psf_fits_name}',f'{psf_fits_path}')

    # GALFIT REQUIREMENTS
    # initial parameters for galfit
    initial_effective_radius,initial_sersic_index,initial_axis_ratio,initial_position_angle = initial_params(galaxy=galaxy_name)

    #initial_param_file_name = f'txt_{galaxy_name}_{version}_initial_params.txt'
    initial_param_file_name = f'txt_{galaxy_name}_initial_params.txt'
    initial_param_file_path = f'{folder_new_version_path}/{initial_param_file_name}'
    initial_param_file = open(f'{initial_param_file_path}','w+')

    # calling the function to create the constraints file
    create_initiaL_params(file=initial_param_file,
            galaxy = galaxy_name,
            version = version,
            eff_rad = initial_effective_radius,
            ser_ind = initial_sersic_index,
            ax_rat = initial_axis_ratio,
            pos_ang = initial_position_angle)

    # closing the constraints file
    initial_param_file.close()


    # constraints for galfit
    n_range,re_range,mag_range,xy_range,q_range,pa_range = constraints_values()

    #constraints_file_name = f'txt_{galaxy_name}_{version}_constraints.txt'
    constraints_file_name = f'txt_{galaxy_name}_constraints.txt'
    constraints_file_path = f'{folder_new_version_path}/{constraints_file_name}'
    constraints_file = open(f'{constraints_file_path}','w+')
    

    # calling the function to create the constraints file
    create_constraints(file=constraints_file,
            n_range = n_range,
            re_range = re_range,
            mag_range = mag_range,
            xy_range = xy_range,
            q_range = q_range,
            pa_range = pa_range)

    # closing the constraints file
    constraints_file.close()

    # OTHER UTILITIES

    if '--version_control' in sys.argv:

        # Create a folder to store the csv files
        csv_folder = f'{galaxy_name}_{version}_csv'
        csv_folder_path = f'{folder_new_version_path}/{csv_folder}'
        if os.path.isdir(f'{csv_folder_path}') == True:
            
            shutil.rmtree(f'{csv_folder_path}')

        os.mkdir(f'{csv_folder_path}')

    # Default version
    if len(sys.argv) == 1:

        # By default the running version is the 
        # medfilt_inout
        sys.argv.append('medfilt_inout')
    
    main()

    print('The anlysis has finished\n')

    # Computing the required time
    end_time = time.time()
    total_time = end_time - start_time
    print(f'The Galfit computing time was {(total_time/3600):1.2f} hours\n')    

    time_file = open(f'{folder_new_version_path}/txt_{galaxy_name}_total_time.txt','w+')
    time_file.write(f'{total_time:.2f} # seconds\n')
    time_file.write(f'{(total_time/60):.2f} # minutes\n')
    time_file.write(f'{(total_time/3600):.2f} # hours\n')

    time_file.close()        
        