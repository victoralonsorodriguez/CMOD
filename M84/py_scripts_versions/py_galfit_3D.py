
# This script contains the median filter alobg with the inout version
# for the analysis of Galfit for a galaxy

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import numpy as np
import pandas as pd
import scipy as scp
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import sys
import os
import subprocess
import shutil
import re

import pdb
import time

from py_galfit_script import create_script
from py_galfit_constraints_script import create_constraints
from py_galfit_constraints_values import constraints_values
from py_galfit_initial_params import initial_params
from py_galfit_initial_params_script import create_initiaL_params
from py_galfit_kpc_correction import kpc_correction


#--------CODE--------#

# to compute the total time
start_time = time.time()

# INFORMATION FROM FOLDER NAME

# getting the current working directory
cwd = os.getcwd()
previus_cwd = cwd + '/..'

# Galaxy's name
galaxy = cwd.split('/')[-2]

# Path to the folder with the datacubes
fits_path = f'{previus_cwd}/{galaxy}_original_fits'

# Checking the latest version run
version_file_path = f'{previus_cwd}/{galaxy}_version.txt'
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
if 'psf' not in folder_name_elements:

    psf_size = 'large'

else:

    psf_size = 'small'


# Name of the psf
psf_fits_name = f'psf_{galaxy}_{psf_size}_flux_norm.fits'
psf_fits_path = f'{folder_new_version_path}/{psf_fits_name}'

# If an specific psf is not created for the galaxy
# it will be used a general small
if os.path.exists(f'{psf_fits_path}') == False:

    psf_fits_name = 'psf_small_stacked_original_flux_norm.fits'
    shutil.copyfile(f'./{psf_fits_name}',f'{psf_fits_path}')


# Obtain fits files names form a directory
fits_list = []
fits_name_list = []

# defining some given data
z_max = 5.4 # maximun redshift
z_step = 0.05 # redshift step interval

# calibrations constants
inst = 0.2 # instrumental scale in arcsec/pix
cal = 25 # calibration constant

# crop factor to find the centre
crop_factor = 10

# GALFIT REQUIREMENTS
# initial parameters for galfit
initial_effective_radius,initial_sersic_index,initial_axis_ratio,initial_position_angle = initial_params()

initial_param_file_name = f'txt_{galaxy}_{version}_initial_params.txt'
initial_param_file_path = f'{folder_new_version_path}/{initial_param_file_name}'
initial_param_file = open(f'{initial_param_file_path}','w+')

# calling the function to create the constraints file
create_initiaL_params(file=initial_param_file,
        galaxy = galaxy,
        version = version,
        eff_rad = initial_effective_radius,
        ser_ind = initial_sersic_index,
        ax_rat = initial_axis_ratio,
        pos_ang = initial_position_angle)

# closing the constraints file
initial_param_file.close()


# constraints for galfit
n_range,re_range,mag_range,xy_range,q_range,pa_range = constraints_values()

constraints_file_name = f'txt_{galaxy}_{version}_constraints.txt'
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

# Create a folder to store the csv files
csv_folder = f'{galaxy}_{version}_csv'
csv_folder_path = f'{folder_new_version_path}/{csv_folder}'
if os.path.isdir(f'{csv_folder_path}') == True:
    
    shutil.rmtree(f'{csv_folder_path}')

os.mkdir(f'{csv_folder_path}')

# Default version
if len(sys.argv) == 1:

    # By default the running version is the 
    # medfilt_inout
    sys.argv.append('medfilt_inout')


#--------STARTING THE ANALYSIS------#

# obtain every needed fits from the directory
for file in sorted(os.listdir(fits_path)):

    # checking for the correct files
    if '.fits' in file and 'temporal' not in file  and 'counts' not in file and 'psf' not in file and 'output' not in file and 'buf' not in file and 'crop' not in file:
    
        fits_list.append(file)
    
        # separation between the name and the extension
        fits_name = file.split('.')[0]
        fits_name_list.append(fits_name)
        
        # create a directory with the same name in order to be keep the file order
        if os.path.isdir(f'{folder_new_version_path}/{fits_name}') == False:
            os.mkdir(f'{folder_new_version_path}/{fits_name}')
        
        shutil.copyfile(f'{fits_path}/{file}',f'{folder_new_version_path}/{file}')

        # open the file in flux
        hdu = fits.open(f'{fits_path}/{file}')
        
        hdr = hdu[0].header # image header for all the subframes
        img = hdu[0].data # image in magnitudes

        # to obtain the flux in counts
        flux = 10**((img - 2.5*np.log10(inst**2)-cal)/-2.5)
        temporal_fits  = f'{fits_name}_temporal.fits'
        temporal_fits_path = f'{folder_new_version_path}/{temporal_fits}'
        
        # saving the fits in counts as temporal intermediate file
        output = fits.writeto(f'{temporal_fits_path}', flux, header=hdr, overwrite=True)
        
        # open the fits in counts
        hdu = fits.open(temporal_fits_path)
    
        img_flux = hdu[0].data # image in counts


        #--------CREATING A DATAFRAME TO STORE THE DATA------#

        # creating the Pandas dataframe to store the information
        df = pd.DataFrame(columns = ['z','X_cent', 'X_cent_err', 'Y_cent','Y_cent_err', 'Mag','Mag_err','Eff_rad','Eff_rad_err','Eff_rad_arcsec','Eff_rad_arcsec_err','kpc_per_arcsec','Eff_rad_kpc','Eff_rad_kpc_err','Ind','Ind_err','Ax_rat','Ax_rat_err','Pos_ang','Pos_ang_err','Chi2','Chi2nu'])
        
        # output path as useful variable
        output_path = f'{folder_new_version_path}/{fits_name}'

        # Export dataframe to text file
        csv_file = open(f'{folder_new_version_path}/{fits_name}_galfit_output.csv','w+')
        


        #------ANALYSING EACH FRAME OF THE DATA CUBE-----#
    
        # loop for all frames layers that form the fits file
        for sub_frame_index in range(0,int(z_max/z_step)+1):
        #for sub_frame_index in range(0,2): # this range just to do a quick 
            
            # Obtaining the corresponding redshift for each image
            redshift_index = sub_frame_index*z_step
            temporal_fits  = f'{fits_name}_z{redshift_index:01.2f}_temporal.fits'
            temporal_fits_path = f'{folder_new_version_path}/{temporal_fits}'

            # obtaining a frame of the fits
            sub_frame = img_flux[sub_frame_index]
            
            # deleting de padding 
            padding = sub_frame == sub_frame[0,0]
            sub_frame[padding] = float('NaN')

            # Depending on the used version
            if sys.argv[1] == 'medfilt_inout' or sys.argv[1] == 'medfilt':
                
                # Median filter applied to each subframe
                sub_frame_med = scp.signal.medfilt(sub_frame, kernel_size = 3)

                # saving the subframe
                output = fits.writeto(f'{temporal_fits_path}', sub_frame_med, header=hdr, overwrite=True)
            
            elif sys.argv[1] == 'inout' or sys.argv[1] == 'original':

                # saving the subframe
                output = fits.writeto(f'{temporal_fits_path}', sub_frame, header=hdr, overwrite=True)

            
            # open the fits in counts without padding
            hdu_galfit = fits.open(f'{temporal_fits_path}')
            img_galfit = hdu_galfit[0].data
            
            # obtaining the shape of the fits
            img_galfit_shape = img_galfit.shape
            nrow = img_galfit.shape[0]
            ncol = img_galfit.shape[1]
            
            nrow_range = (nrow//crop_factor)
            ncol_range = (ncol//crop_factor)
            
            # obtaining the nearest to the center maximun pixel value position to center the model
            img_galfit_center = img_galfit[nrow//2 - nrow_range:nrow//2 + nrow_range, ncol//2 - ncol_range:ncol//2 + ncol_range]
            
            #fits.writeto(f'{sub_frame_index}_crop.fits', img_galfit_center, header=hdr, overwrite=True)
            
            max_pix_value_center = np.nanmax(img_galfit_center)
            max_pix_value_pos_x = np.where(img_galfit == max_pix_value_center)[0][0]+1
            max_pix_value_pos_y = np.where(img_galfit == max_pix_value_center)[1][0]+1
            
            max_pix_pos = (max_pix_value_pos_x ,max_pix_value_pos_y)
            
            # obtaining the surface magnitude of the image
            total_flux = np.nansum(img_galfit)
            img_magnitude = round(-2.5*np.log10(total_flux)+cal,2)
            
            # creating the script
            redshift_index = sub_frame_index*z_step
            script = f'{fits_name}_z{redshift_index:01.2f}_temporal.script'
            script_path = f'{folder_new_version_path}/{script}'
            f = open(f'{script_path}','w+')
            
            # the ouput file name will be
            output_file = fits_name +f'_z{redshift_index:01.2f}' +'_galfit_output.fits'
            output_file_path = f'{folder_new_version_path}/{output_file}'

            # closing the constraints file
            constraints_file.close()

            # For versions using inout methodology
            if sys.argv[1] == 'inout' or sys.argv[1] == 'medfilt_inout':
    
                # For all frames except for the first one we introduce the output 
                # parameters of previous frame as input
                if sub_frame_index != 0:

                    # calling the function to create the script
                    create_script(file=f,input_file_name=temporal_fits,
                            output_file_name = output_file,
                            file_shape = img_galfit_shape,
                            zp_const = cal,
                            pix_scl = inst,
                            img_center = (y_center,x_center),
                            int_mag = mag,
                            eff_rad = eff_rad,
                            ser_ind = ser_index,
                            ax_rat = ax_rat,
                            pos_ang = pos_ang,
                            psf= psf_fits_name,
                            cons_file = constraints_file_name)

                # Initial parameters for the first frame of the filter datacube
                elif sys.argv[1] == 'medfilt' or sys.argv[1] == 'original':

                    # calling the function to create the script
                    create_script(file=f,input_file_name=temporal_fits,
                            output_file_name = output_file,
                            file_shape = img_galfit_shape,
                            zp_const = cal,
                            pix_scl = inst,
                            img_center = max_pix_pos,
                            int_mag = img_magnitude,
                            eff_rad = initial_effective_radius,
                            ser_ind = initial_sersic_index,
                            ax_rat = initial_axis_ratio,
                            pos_ang = initial_position_angle,
                            psf= psf_fits_name,
                            cons_file = constraints_file_name)
                    
            
            # For versions not using inout methodology
            elif sys.argv[1] == 'original' or sys.argv[1] == 'medfilt':

                # Initial parameters for every frame of the filter datacube
                create_script(file=f,input_file_name=temporal_fits_path,
                        output_file_name = output_file_path,
                        file_shape = img_galfit_shape,
                        zp_const = cal,
                        pix_scl = inst,
                        img_center = max_pix_pos,
                        int_mag = img_magnitude,
                        eff_rad = initial_effective_radius,
                        ser_ind = initial_sersic_index,
                        ax_rat = initial_axis_ratio,
                        pos_ang = initial_position_angle,
                        psf= psf_fits_name,
                        cons_file = constraints_file_name)
    

            # closing the script file
            f.close()

            # We need to change the cwd because galfit save the fit.log
            # and the other data in the cwd
            os.chdir(f'{folder_new_version_path}')

            # runing galfit for the subframe fits with the corresponding
            subprocess.run(['galfit', script_path])

            # Changin again to the original cwd
            os.chdir(f'{cwd}')

            #-----OBTAINING THE PARAMETERS OF THE ANALYSIS----#
            
            # Extracting the galfit results data from the header of each output model frame
            
            # This 'if' is here because sometimes galfit doesn't fit a parameter and due that
            # it doesn't create an output file. This 'if' is to avoid the error of file not found
            if os.path.isfile(f'{output_file_path}'):
                
                # open the fits
                hdu = fits.open(f'{output_file_path}')
                
                # the output model is the layer 2
                hdr = hdu[2].header 
                img = hdu[2].data
                
                # find the redshif in the ouput file name
                z = re.findall('z\d.\d\d',output_file_path)[0][1:]

                #pdb.set_trace()
                
                # output values from galfit and their errors
                x_center = hdr['1_XC'].split(' ')[0]
                x_center_err = hdr['1_XC'].split(' ')[2]
                
                y_center = hdr['1_YC'].split(' ')[0]
                y_center_err = hdr['1_YC'].split(' ')[2]
                
                mag = hdr['1_MAG'].split(' ')[0]
                mag_err = hdr['1_MAG'].split(' ')[2]

                if '*' in mag:
                    mag = np.nan
                    mag_err = np.nan
                
                eff_rad = hdr['1_RE'].split(' ')[0]
                eff_rad_err = hdr['1_RE'].split(' ')[2]

                if '*' in eff_rad:
                    eff_rad = np.nan
                    eff_rad_err = np.nan
                
                # Effective radius in arcsec
                eff_rad_arcsec = float(eff_rad) * inst
                eff_rad_arcsec_err = float(eff_rad_err) * inst
                
                # Effective radius in kpc
                #moving distances from kpc to arcsec
                cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
                kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(z)
                
                # Extract the unit to create a Quantity inside the if
                kpc_unit = kpc_per_arcsec.unit

                '''
                eff_rad_kpc = eff_rad_arcsec * kpc_per_arcsec.value    
                eff_rad_kpc_err = eff_rad_arcsec_err * kpc_per_arcsec.value
                
                if z == '0.00':
                    
                    # Obtaining a Quantity
                    kpc_per_arcsec = kpc_correction(galaxy)*kpc_unit
                    
                    eff_rad_kpc = eff_rad_arcsec * kpc_per_arcsec.value              
                    eff_rad_kpc_err = eff_rad_arcsec_err * kpc_per_arcsec.value 

                '''
                
                kpc_per_arcsec = kpc_correction(galaxy)*kpc_unit
                    
                eff_rad_kpc = eff_rad_arcsec * kpc_per_arcsec.value              
                eff_rad_kpc_err = eff_rad_arcsec_err * kpc_per_arcsec.value 
                
                ser_index = hdr['1_N'].split(' ')[0]
                ser_index_err = hdr['1_N'].split(' ')[2]


                if '*' in ser_index:
                    ser_index = np.nan
                    ser_index_err = np.nan
                
                ax_rat = hdr['1_AR'].split(' ')[0]
                ax_rat_err = hdr['1_AR'].split(' ')[2]

                if '*' in ax_rat:
                    ax_rat = np.nan
                    ax_rat_err = np.nan
                
                pos_ang = hdr['1_PA'].split(' ')[0]
                pos_ang_err = hdr['1_PA'].split(' ')[2]

                if '*' in pos_ang:
                    pos_ang = np.nan
                    pos_ang_err = np.nan
                
                chi2 = hdr['CHISQ']
                chi2nu = hdr['CHI2NU']
                
                # adding all the values to a new column of the dataframe
                df.loc[len(df)]=[z,x_center,x_center_err,y_center,y_center_err,mag,mag_err,eff_rad,eff_rad_err,eff_rad_arcsec,eff_rad_arcsec_err,kpc_per_arcsec.value,eff_rad_kpc,eff_rad_kpc_err,ser_index,ser_index_err,ax_rat,ax_rat_err,pos_ang,pos_ang_err,chi2,chi2nu]
                df = df.astype('float64')


        #----WRITING THE DATAFRAME INTO A CSV----#

        # Sort the dataframe by redshift
        df.sort_values(by=['z'],ascending=True, inplace=True, ignore_index=True)
        
        df_string = df.to_csv(header=True, index=False, sep=',')

        csv_file.write(df_string)

        csv_file.close()


        #----ORGANIZING ALL THE CREATED FILES INTO NEW FOLRDES----#
        
        # moving all the output and created files to a folder in order to keep the order in the directory        
        for output_file in sorted(os.listdir(folder_new_version_path)):

            if '.csv' in output_file:

                shutil.copyfile(f'{folder_new_version_path}/{output_file}',f'{folder_new_version_path}/{csv_folder}/{output_file}')

            if '.py' not in output_file and ('galfit.' in output_file or '_galfit_output' in output_file or 'fit.' in output_file or 'temporal' in output_file):
                
                os.replace(f'{folder_new_version_path}/{output_file}', f'{folder_new_version_path}/{fits_name}/{output_file}')

        
        # copy the original fits file to the new folder
        os.remove(f'{folder_new_version_path}/{file}')
        


print('The anlysis has finished\n')

# Computing the required time
end_time = time.time()
total_time = end_time - start_time
print(f'The Galfit computing time was {(total_time/3600):1.2f} hours\n')    

time_file = open(f'{folder_new_version_path}/txt_{galaxy}_total_time.txt','w+')
time_file.write(f'{total_time:.2f} # seconds\n')
time_file.write(f'{(total_time/60):.2f} # minutes\n')
time_file.write(f'{(total_time/3600):.2f} # hours\n')

time_file.close()        
        

        
        
        
        
        
        
        
        
        
