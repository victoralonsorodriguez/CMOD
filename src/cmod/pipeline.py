import numpy as np
import pandas as pd
import scipy as scp
from pathlib import Path

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import os
import re
import sys
import subprocess
import shutil
import mmap
import struct
import traceback

import pdb
import time
import warnings
warnings.filterwarnings('ignore')

from .cosmology import zlocal_correction, z_counts_dimming, z_dimming_factor
from .fitting import (constraints_values, create_constraints, 
                          create_script, create_conf_imfit,
                          initial_params, create_initiaL_params)
from .io import open_fits, argparse_values, extract_filter_name
from .photometry import fits_mag_to_counts, values_counts_to_mag, filter_pixel_scale_dict
from .processing import max_center_value, isophote_fitting, initial_conditions, create_psf
from .resampling import resampling_frame_pixels,resampling_frame, zoom_factor
from .results import galfit_init_dataframe, galfit_create_dataframe,imfit_init_dataframe,imfit_create_dataframe
from .utils import (Cronometro, round_number, create_folder, rad_to_deg_abs,
                        version_directory, version_file, version_file_last)


def run_analysis_pipeline():
    
    #############################################################
    ###------------------GENERAL INFORMATION------------------###
    #############################################################
    
    cronometro = Cronometro()
    cronometro.iniciar()
    
    # Computing the total analysis time
    start_time = time.time()
    
    # Loading configuration
    (analysis_programme,analysis_version,analysis_range,
                version_continuation,version_control,psf,initial_conditions_mode,
                ori_pixel_scale,zcal_const,crop_factor) = argparse_values(phase='analysis')
    
    # Obtaining the current working directory
    # in which the programme is executed
    PATH_FILE = Path(__file__).resolve()
    PATH_SRC_CMOD = PATH_FILE.parent
    PATH_SRC = PATH_SRC_CMOD.parent
    PATH_REPO = PATH_SRC.parent
    PATH_DATA = PATH_REPO / "data"
    PATH_OUTPUTS = PATH_REPO / "outputs"
    PATH_LOGS = PATH_OUTPUTS / "logs"
    
    datacube_folder_name_list = []
    datacube_folder_path_list = []
    
    for path, subdir, file in sorted(os.walk(PATH_DATA)):
        for subdir_name in subdir:
            if 'original_fits' in subdir_name:
                datacube_folder_name_list.append(subdir_name)
                datacube_folder_path_list.append(path)
    
    for path_pos,path in enumerate(datacube_folder_path_list): 
        
        # Changing the directory to obtain the original path
        # os.chdir(path)
        # cwd_galaxy = os.getcwd()
        
        # Returning to the scripts directory
        # os.chdir(PATH_SRC_CMOD)

        # Galaxy name obtaiane from datacube directories
        galaxy_name = (datacube_folder_name_list[path_pos].split('/')[-1]).split('_')[0]
        
        # Absolute path to the original fits datacubes
        # This files are the ones to be analyzed
        datacube_original_folder_path = f'{PATH_DATA}/{galaxy_name}_original_fits'
        
        # Obtaining the datacube names
        datacube_name_list = []
        for datacube_name in sorted(os.listdir(datacube_original_folder_path)):
            if '.fits' in datacube_name:
                datacube_name_list.append(datacube_name)
        
        # Datacubes general name for the galaxy without the filter name
        datacube_name_general = '_'.join((datacube_name_list[0].split('.fits')[0]).split('_')[:-1])
        
        # We need the  pixel scale and the calibration constant to convert 
        # from counts to magnitudes
        pxsc_zcal_const = (ori_pixel_scale,zcal_const)
        
        
        ###------------------VERSION CONTROL------------------###
        #########################################################
        
        # Galaxy analysis new directory
        next_version = '00'
        version_base_name = f'{galaxy_name}_analysis'
        version_file_name = f'txt_{galaxy_name}_version.txt'
        version_file_path = f'{PATH_DATA}/{version_file_name}'
        
        # If version control is activated
        if version_control == 1:
            
            # Obtained the last runned version from the
            # version file and the directory version name
            prev_ver_file = version_file_last(version_file_path)
            prev_ver_dir = version_directory(cwd=PATH_DATA,
                                            base_name = version_base_name,
                                            original_dir = PATH_DATA)
            
            next_version = max(prev_ver_file,prev_ver_dir)+1
            if next_version < 10:
                next_version = f'0{next_version}'
                
        # Creating the version file and the version directory
        version_file(galaxy_name,next_version,analysis_version,psf,
                    archivo=version_file_path,
                    analysis_programe=analysis_programme)
        analysis_folder_name = f'{version_base_name}_V{next_version}_{analysis_programme}'
        analysis_folder_path = f'{PATH_DATA}/{analysis_folder_name}'
        create_folder(analysis_folder_path)

        ###------------------POINT SPREAD FUNCTION------------------###
        ###############################################################

        # Name of the psf
        if psf == 'large':    
            
            fwhm = 2.0336 # FWHM of the PSF
            beta = 1.8781 # Beta parameter of the PSF
            
            (psf_fits_name,
             psf_fits_path) =  create_psf(galaxy_name=galaxy_name,
                                            datacube_folder_path=datacube_original_folder_path,
                                            analysis_folder_path=analysis_folder_path,
                                            pixel_scale=ori_pixel_scale,
                                            zcal_const=zcal_const,
                                            fwhm=fwhm,
                                            beta=beta)
            
            os.chdir(PATH_SRC_CMOD)
            
        elif psf == 'small':
            psf_fits_name = f'psf_small_stacked_original_flux_norm.fits'
            psf_fits_path_ori = f'{PATH_REPO}/data/{psf_fits_name}'
            psf_fits_path = f'{analysis_folder_path}/{psf_fits_name}'
            shutil.copyfile(f'{psf_fits_path_ori}',f'{psf_fits_path}')
        else:
            psf_fits_name = psf.split('/')[-1]
            psf_fits_path = psf
            

        ###------------------INITIAL CONDITIONS------------------###
        ############################################################

        # GALFIT REQUIREMENTS
        # initial parameters for galfit
        initial_effective_radius,initial_sersic_index,initial_axis_ratio,initial_position_angle = initial_params(galaxy=galaxy_name)

        #initial_param_file_name = f'txt_{galaxy_name}_{version}_initial_params.txt'
        initial_param_file_name = f'txt_{galaxy_name}_initial_params.txt'
        initial_param_file_path = f'{analysis_folder_path}/{initial_param_file_name}'
        initial_param_file = open(f'{initial_param_file_path}','w+')
        
        # calling the function to create the constraints file
        create_initiaL_params(file=initial_param_file,
                galaxy = galaxy_name,
                version = analysis_version,
                eff_rad = initial_effective_radius,
                ser_ind = initial_sersic_index,
                ax_rat = initial_axis_ratio,
                pos_ang = initial_position_angle)

        # closing the constraints file
        initial_param_file.close()


        ###------------------CONSTRAINTS------------------###
        #####################################################

        # constraints for galfit
        n_range,re_range,mag_range,xy_range,q_range,pa_range = constraints_values()

        #constraints_file_name = f'txt_{galaxy_name}_{version}_constraints.txt'
        constraints_file_name = f'txt_{galaxy_name}_constraints.txt'
        constraints_file_path = f'{analysis_folder_path}/{constraints_file_name}'
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


        ###------------------CSV FOLDER------------------###
        ####################################################    
        
        # CSV folder name
        csv_folder_name = f'{analysis_folder_name}_csv'

        # Create a folder to store the csv files
        # in the new version directory
        csv_folder_path = f'{analysis_folder_path}/{csv_folder_name}'
        create_folder(csv_folder_path,overwrite=False)

        ################################################################
        ###------------------ANALYSYS CODE FUNCTION------------------###
        ################################################################
        
        try:
            print(f'\nAnalysing {galaxy_name} with {analysis_programme} '
                    f'in the run {next_version}')
            analisys(datacube_name_list, analysis_folder_path, analysis_programme, 
                     datacube_original_folder_path, pxsc_zcal_const, analysis_range, 
                     crop_factor, initial_conditions_mode, galaxy_name, 
                     psf_fits_name, constraints_file_name, csv_folder_path,
                     analysis_folder_name,analysis_version,zcal_const,ori_pixel_scale,
                     PATH_SRC_CMOD)
            print(f'\nThe analysis for {galaxy_name} with {analysis_programme} '
                    f'in the run {next_version} has finished')
            
        except Exception as e:
            print(f'\nError in the analysis module:\n{e}')

        
        ###------------------TIME INFORMATION------------------###
        ##########################################################

        # Computing the required time
        
        end_time = time.time()
        total_time = end_time - start_time
        print(f'\n{analysis_programme} computing time was {(total_time/3600):1.2f} hours\n')    

        time_file = open(f'{analysis_folder_path}/txt_{galaxy_name}_total_time.txt','w+')
        time_file.write(f'{total_time:.2f} # seconds\n')
        time_file.write(f'{(total_time/60):.2f} # minutes\n')
        time_file.write(f'{(total_time/3600):.2f} # hours\n')

        time_file.close()        



def analisys(datacube_name_list, analysis_folder_path, analysis_programme, 
             datacube_original_folder_path, pxsc_zcal_const, analysis_range, 
             crop_factor, initial_conditions_mode, galaxy_name, 
             psf_fits_name, constraints_file_name, csv_folder_path,
             analysis_folder_name,analysis_version,zcal_const,ori_pixel_scale,
             PATH_SRC_CMOD):
    
    '''
    This function contains the analysis process
    '''
    
    # The analysis is carried out for each datacube
    for datacube_pos,datacube_name in enumerate(datacube_name_list):
        
        ###------------------MANIPULATING DATACUBE------------------###
        ###############################################################
        
        datacube_name_noext = datacube_name.split('.fits')[0]
        
        # Creating a directory with the datacube name to store the 
        # files produced during the analysis
        datacube_folder_path = f'{analysis_folder_path}/{datacube_name_noext}'
        create_folder(datacube_folder_path,
                      overwrite=False)
        
        # Creating a subfolder inside each datacube folder
        if analysis_programme == 'Galfit':
            datacube_galfit_folder_name = f'{datacube_name_noext}_galfit_output'
            datacube_galfit_folder_path = f'{datacube_folder_path}/{datacube_galfit_folder_name}'
            create_folder(datacube_galfit_folder_path)
        
        # Copying the datacube from the original folder
        # to its corresponding folder
        datacube_original_path = f'{datacube_original_folder_path}/{datacube_name}'
        shutil.copyfile(datacube_original_path,
                        f'{datacube_folder_path}/{datacube_name}')

        # The datacube is in magnitudes but for the analysis is
        # required its value in counts       
        # The images have a padding that should be removed 
        if analysis_programme == 'Galfit':
            pad_val = np.nan
        elif analysis_programme == 'Imfit':
            pad_val = 0
        
        datacube_counts_path =  fits_mag_to_counts(fits_input_path = datacube_original_path,
                                                   fits_output_path = datacube_folder_path,
                                                   const = pxsc_zcal_const,
                                                   remove_padding = True,
                                                   padding_value=pad_val)
        
        # We need to open the fits to analyse each frame
        datacube_flux_hdr, datacube_flux_data, datacube_flux_name  = open_fits(datacube_counts_path)
        
        # Redshift step of the datacube
        redshift_step = datacube_flux_hdr['CDELT3']
        
        # Filter system
        filter_name = extract_filter_name(datacube_name)
        sim_pixel_scale = filter_pixel_scale_dict(filter_name)
        
        
        ###------------------CREATING A DATAFRAME------------------###
        ##############################################################

        # Creating a dataframe to storing analysis data results
        if analysis_programme == 'Galfit':
            df = galfit_init_dataframe()
        elif analysis_programme == 'Imfit':
            df = imfit_init_dataframe()

        # Export dataframe to text file for each datacube
        csv_file_name = f'{datacube_name_noext}_{analysis_programme}.csv'
        csv_file_path = open(f'{analysis_folder_path}/{csv_file_name}','w+')
        
        
        ###------------------FRAME ANALYSIS------------------###
        ########################################################
        
        # Determining frame analysing range
        if analysis_range == 0:
            frame_range = datacube_flux_data.shape[0]
        else:
            frame_range = analysis_range
            
        try:
        
            # Analysing each frame individually
            for frame_pos in range(frame_range):
                
                
                # Obtaining the corresponding redshift of each frame
                redshift_value = round_number(frame_pos * redshift_step,2)
                zlocal_value = zlocal_correction(galaxy=galaxy_name)
                
                # We are going to export each frame of the data cube individually
                frame_ori = datacube_flux_data[frame_pos]
                
                frame_ori_name = f'{datacube_name.split(".")[0]}_z{redshift_value:.2f}_counts.fits'
                frame_ori_name_noext = frame_ori_name.split('.fits')[0]
                frame_ori_path = f'{analysis_folder_path}/{frame_ori_name}'
                
                print(f'\n\n\tAnalyzing {frame_ori_name} with {analysis_programme}\n')     

                # Depending on the used version we are going to use a median 
                # filter in order to reduce the image noise
                # By default the filter is aplied
                if ('--medfilt_inout' in sys.argv or '--medfilt' in sys.argv or
                    '--inout' not in sys.argv or '--original' not in sys.argv):
                    
                    # Median filter applied to each frame                    
                    frame_med = scp.ndimage.median_filter(frame_ori, size=3, cval=np.nan)
                    frame_ori = frame_med

                # saving the subframe
                fits.writeto(frame_ori_path, frame_ori, header=datacube_flux_hdr, overwrite=True)
                
                # Obtaining the frame shape
                frame_ori_shape = frame_ori.shape
                
                # Resampling the frame
                frame_res_shape = resampling_frame_pixels(galaxy=galaxy_name,
                                                               galaxy_shape = frame_ori_shape,
                                                               original_sensor_pixscale=ori_pixel_scale,
                                                               zsim=redshift_value,
                                                               simulated_sensor_pxscale=sim_pixel_scale)
                
                frame_res_path = resampling_frame(frame_path=frame_ori_path,
                                                  resample_frame_size=frame_res_shape)
                
                
                
                # Opening resampling fits
                hdr_res,frame_res,frame_res_name_noext = open_fits(frame_res_path)
                frame_res_name = f'{frame_res_name_noext}.fits'
                
                zoom_applied = zoom_factor(frame_ori_shape,frame_res_shape)
                
                if redshift_value != 0.0:
                    frame_dim = frame_res * z_dimming_factor(redshift_value)
                else:
                    frame_dim = frame_res * z_dimming_factor(zlocal_value)
                
                frame_dim_name = f'{frame_res_name_noext}_dim.fits'
                frame_dim_name_noext = frame_dim_name.split('.fits')[0]
                frame_dim_path = f'{analysis_folder_path}/{frame_dim_name}'
                
                frame_dim_shape = frame_res_shape
                
                fits.writeto(frame_dim_path, frame_dim, header=hdr_res, overwrite=True)
                
                # Obtaining the center of the galaxy as the maximun central pixel value
                # Only for the first image
                frame_dim_max_pix_value, frame_dim_max_pix_value_pos = max_center_value(frame_dim,crop_factor)

                #if frame_pos == 0:
                    
                try:
                    if initial_conditions_mode == 0:
                        gal_iso_fit_csv_path = isophote_fitting(frame_dim_path,
                                                                frame_dim_max_pix_value_pos,
                                                                frame_dim_shape,
                                                                cons=(sim_pixel_scale,zcal_const),
                                                                output_path=analysis_folder_path)
                        gal_iso_fit_df = pd.read_csv(gal_iso_fit_csv_path)
                        (break_pos, break_pos_bul, I_0_disk_in, 
                            h_disk_in, n_bul_in, r_e_bul_in, I_e_bul_in) = initial_conditions(frame_dim_name_noext,
                                                                                            gal_iso_fit_df, 
                                                                                            'sma', 
                                                                                            'intens', 
                                                                                            'intens_err',
                                                                                            const=(sim_pixel_scale,zcal_const),
                                                                                            output_path=analysis_folder_path)
                            
                        # Computing the mean values from the external partes of the galaxy
                        ell_mean_in = round_number(gal_iso_fit_df.loc[break_pos:, 'ellipticity'].mean(),3)
                        ax_rat_in = round_number(1 - ell_mean_in,3)
                        pa_rad_mean_in = gal_iso_fit_df.loc[break_pos:, 'pa'].mean()
                        pa_deg_in = round_number(rad_to_deg_abs(pa_rad_mean_in)-90,3)
                        
                    else:
                        r_e_bul_in,n_bul_in,ax_rat_in,pa_deg_in = initial_params(galaxy=galaxy_name)

                except Exception as e:
                    print(f'Following error obtaining initial conditions'
                        f'\n{e}')

                
                # Obtaining them toal surface magnitude of the frame
                frame_dim_counts = np.nansum(frame_dim)
                
                frame_dim_total_mag = values_counts_to_mag(frame_dim_counts,
                                                            const = pxsc_zcal_const)    
                
                # Creating a script to run galfit or imfit
                if analysis_programme == 'Galfit':
                    script_name = f'{frame_res_name_noext}_{analysis_programme}.script'
                elif analysis_programme == 'Imfit':
                    script_name = f'{frame_res_name_noext}_{analysis_programme}.txt'
                    
                script_path = f'{analysis_folder_path}/{script_name}'
                script_file = open(f'{script_path}','w+')

                
                # the ouput file name will be
                output_file_name = f'{frame_dim_name_noext}_{analysis_programme}.fits'
                output_file_path = f'{analysis_folder_path}/{output_file_name}'
            
                residual_name = f'{frame_dim_name_noext}_{analysis_programme}_residual.fits'
                residual_path = f'{analysis_folder_name}/{residual_name}'
                
                # For versions using inout methodology
                if analysis_version == 'inout' or analysis_version == 'medinout':
        
                    # For all frames except for the first one we introduce the output 
                    # parameters of previous frame as input
                    if frame_pos != 0:
                        if analysis_programme == 'Galfit':
                            
                            # calling the function to create the script
                            create_script(file=script_file,
                                        input_file_name=frame_res_name,
                                        output_file_name = output_file_name,
                                        file_shape = frame_res_shape,
                                        zp_const = zcal_const,
                                        pix_scl = sim_pixel_scale,
                                        img_center = frame_dim_max_pix_value_pos,
                                        int_mag = frame_dim_total_mag,
                                        eff_rad = r_e_bul_in,
                                        ser_ind = n_bul_in,
                                        ax_rat = ax_rat_in,
                                        pos_ang = pa_deg_in,
                                        psf= psf_fits_name,
                                        cons_file = constraints_file_name)
                            
                        elif analysis_programme == 'Imfit':
                            if I_e <= 1:
                                I_e += 20
                            if pos_ang <= 1:
                                pos_ang = 2 
                            if eff_rad <= 1:
                                eff_rad = 2
                                
                            create_conf_imfit(file=script_file,
                                    funct=['Sersic'],
                                    galaxy = galaxy_name,
                                    img_center_x = [x_center,x_center-1,x_center+1],
                                    img_center_y = [y_center,y_center-1,y_center+1],
                                    pa_ser = [pos_ang,1,179],
                                    ell_ser= [ell,0.05,0.95],
                                    n=[ser_index,0.6,6],
                                    r_e=[eff_rad,1,200],
                                    I_e=[I_e,1,
                                        round_number(I_e+50,3)])

                    # Initial parameters for the first frame of the filter datacube
                    else:
                        if analysis_programme == 'Galfit':
                            # calling the function to create the script
                            create_script(file=script_file,
                                        input_file_name=frame_res_name,
                                        output_file_name = output_file_name,
                                        file_shape = frame_res_shape,
                                        zp_const = zcal_const,
                                        pix_scl = sim_pixel_scale,
                                        img_center = frame_dim_max_pix_value_pos,
                                        int_mag = frame_dim_total_mag,
                                        eff_rad = r_e_bul_in,
                                        ser_ind = n_bul_in,
                                        ax_rat = ax_rat_in,
                                        pos_ang = pa_deg_in,
                                        psf= psf_fits_name,
                                        cons_file = constraints_file_name)
                        
                        elif analysis_programme == 'Imfit':
                            if pa_deg_in < 0:
                                pa_deg_in += 90
                            elif pa_deg_in >180:
                                pa_deg_in -= 90
                            pa_deg_in = round_number(pa_deg_in,3)
                            
                            if I_e_bul_in <= 1:
                                I_e_bul_in = round_number(frame_dim_max_pix_value/5, 3)

                            
                            create_conf_imfit(file=script_file,
                                    funct=['Sersic'],
                                    galaxy = galaxy_name,
                                    img_center_x = [frame_dim_max_pix_value_pos[1],frame_dim_max_pix_value_pos[1]-1,frame_dim_max_pix_value_pos[1]+1],
                                    img_center_y = [frame_dim_max_pix_value_pos[0],frame_dim_max_pix_value_pos[0]-1,frame_dim_max_pix_value_pos[0]+1],
                                    pa_ser = [pa_deg_in,0,180],
                                    ell_ser= [ell_mean_in,0.05,0.95],
                                    n=[n_bul_in,0.6,6],
                                    r_e=[r_e_bul_in,0,200],
                                    I_e=[round_number(I_e_bul_in,3),
                                        1,round_number(I_e_bul_in+50,3)])

                # For versions not using inout methodology
                elif 'original' in sys.argv or 'medfilt' in sys.argv:
                    if analysis_programme == 'Galfit':
                        # Initial parameters for every frame of the filter datacube
                        create_script(file=script_file,
                                    input_file_name=frame_res_name,
                                    output_file_name = output_file_name,
                                    file_shape = frame_res_shape,
                                    zp_const = zcal_const,
                                    pix_scl = sim_pixel_scale,
                                    img_center = frame_dim_max_pix_value_pos,
                                    int_mag = frame_dim_total_mag,
                                    eff_rad = r_e_bul_in,
                                    ser_ind = n_bul_in,
                                    ax_rat = ax_rat_in,
                                    pos_ang = pa_deg_in,
                                    psf= psf_fits_name,
                                    cons_file = constraints_file_name)
                        
                    elif analysis_programme == 'Imfit':
                            pass

                
                # closing the script file
                script_file.close()

                os.chdir(analysis_folder_path)
                
                if analysis_programme == 'Galfit':

                    # runing galfit for the subframe fits with the corresponding
                    subprocess.run(['galfit', script_name])
                    
                elif analysis_programme == 'Imfit':
                    
                    subprocess.run(['imfit',
                                    frame_res_name,
                                    '-c', script_name,
                                    #'--mask',mask_path,
                                    '--gain=','4.76',
                                    '--readnoise=','5.0625',
                                    '--sky','100',
                                    '--psf',psf_fits_name,
                                    '--save-model=',output_file_name,
                                    '--save-residual=',residual_name])     
                
                
                
                if analysis_programme == 'Imfit':
                    
                    fitting_file_name_ori = 'bestfit_parameters_imfit.dat'
                    fitting_file_name_new = f'{frame_res_name_noext}_{analysis_programme}_bestfit_param.dat'
                    fitting_file_path_new = f'{analysis_folder_path}/{fitting_file_name_new}'
                    if os.path.isfile(f'{fitting_file_name_ori}'):
                        shutil.move(fitting_file_name_ori, fitting_file_name_new)
                    
                os.chdir(PATH_SRC_CMOD)  
                
            
                #-----OBTAINING THE PARAMETERS OF THE ANALYSIS----#
                
                # Extracting the galfit results data from the header of each output model frame
                
                # This 'if' is here because sometimes galfit doesn't fit a parameter and due that
                # it doesn't create an output file. This 'if' is to avoid the error of file not found
                if analysis_programme == 'Galfit':
                
                    if os.path.isfile(f'{output_file_path}'):
                        
                        df,y_center,x_center,mag,eff_rad,ser_index,ax_rat,pos_ang = galfit_create_dataframe(df,galaxy_name,output_file_path,pxsc_zcal_const)

                        
                elif analysis_programme == 'Imfit':
                    if os.path.isfile(f'{fitting_file_path_new}'):
                        
                        df,y_center,x_center,I_e,eff_rad,ser_index,ell,pos_ang = imfit_create_dataframe(df=df,
                                                                                                        galaxy_name=galaxy_name,
                                                                                                        redshift=redshift_value,
                                                                                                        param_file=fitting_file_path_new,
                                                                                                        pxsc_zcal_const=pxsc_zcal_const)
                    
                    
                
                    
            #----WRITING THE DATAFRAME INTO A CSV----#

            # Sort the dataframe by redshift
            df.sort_values(by=['z'],ascending=True, inplace=True, ignore_index=True)
            
            df_string = df.to_csv(header=True, index=False, sep=',')

            csv_file_path.write(df_string)

            csv_file_path.close()
            
        except Exception as e:
            
            error_traceback = traceback.format_exc()
            
            print(f"\n\n--- ERROR during filter analysis ---")
            print(f"Error Type: {type(e).__name__}") 
            print(f"Error Message: {e}")           
            print(f"Traceback:\n{error_traceback}")
            
            
        #----ORGANIZING ALL THE CREATED FILES INTO NEW FOLRDES----#
        
        # moving all the output and created files to a folder in order to keep the order in the directory        
        for output_file in sorted(os.listdir(analysis_folder_path)):
            
            file_original_path = f'{analysis_folder_path}/{output_file}'
            
            if '.csv' in output_file and analysis_programme in output_file:
                file_destiny_path = f'{csv_folder_path}/{output_file}'
                
                shutil.copyfile(file_original_path,file_destiny_path)

            if ('.py' not in output_file and os.path.isdir(output_file) == False and 
                (datacube_name_noext in output_file or 'galfit.' in output_file or 'fit.log' in output_file)):
                
                if analysis_programme == 'Galfit':  
                    
                    if 'galfit.' in output_file or 'fit.log' in output_file:

                        file_destiny_path = f'{datacube_galfit_folder_path}/{output_file}'
                        os.replace(file_original_path, file_destiny_path)
                        continue
                
                file_destiny_path = f'{datacube_folder_path}/{output_file}'
                if file_original_path != '/'.join(file_destiny_path.split('/')[:-1]):
                    
                    os.replace(file_original_path, file_destiny_path)