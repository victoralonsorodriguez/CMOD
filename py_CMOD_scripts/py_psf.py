
# This script creates a PSF of the size of the galaxy image

from astropy.io import fits
import numpy as np

import os
import subprocess

import pdb

from py_psf_script import create_psf_script
from cmod.io import open_fits, argparse_values


'''
Hay que cambiar este código para que tenga en cuenta los flags en la linea de comandos 
y además para que la PSF sea un pixel más grande que la imagen en caso de que la imagen
de la galaxia no sea impar para que quede completamente centrada
'''


def create_psf(galaxy_name,
               datacube_folder_path,
               analysis_folder_path,
               pixel_scale,
               zcal_const,
               fwhm,
               beta):
    
    datacube_name_list = []
    
    for datacube in sorted(os.listdir(datacube_folder_path)):
        if '.fits' in datacube:
            datacube_name_list.append(datacube)
            

    datacube_psf_name = datacube_name_list[0]
    datacube_psf_path = f'{datacube_folder_path}/{datacube_psf_name}'

    _,datacube_img,_ = open_fits(datacube_psf_path)
    
    frame_psf_img = datacube_img[0]
    
    # obtaining the shape
    x_len = frame_psf_img.shape[0]
    y_len = frame_psf_img.shape[1]

    
    if x_len % 2 == 0:
        x_len += 1
    if y_len % 2 == 0:
        y_len += 1
    
    psf_shape = (y_len,x_len)
    psf_center = (y_len//2 + 1,x_len//2 + 1)
    
    # creating the script for running galfit
    psf_script_name = f'psf_{galaxy_name}_large.script'
    psf_script_path = f'{analysis_folder_path}/{psf_script_name}'
    
    psf_script_file = open(psf_script_path,'w+')
    
    psf_output_name = f'psf_{galaxy_name}_large.fits'
    psf_output_path = f'{analysis_folder_path}/{psf_output_name}'
    
    create_psf_script(file=psf_script_file,
                input_file_name='none',
                output_file_name = psf_output_path,
                psf = 'none',
                cons_file = 'none',
                file_shape = psf_shape,
                zp_const = zcal_const,
                pix_scl = pixel_scale,
                img_center = psf_center,
                int_mag = 1,
                fwhm = fwhm,
                beta = beta,
                ax_rat = 1,
                pos_ang = 0)
    
    psf_script_file.close()
    
    # creating the PSF with the given parameters
    os.chdir(analysis_folder_path)
    subprocess.run(['galfit', psf_script_name])
    
    psf_hdr,psf_img,_ = open_fits(psf_output_path)
    
    # total flux of the PSF
    total_flux = np.sum(psf_img)
    
    # normalizing the PSF
    psf_norm = psf_img / total_flux

    # saving the result
    psf_norm_name = f'psf_{galaxy_name}_large_flux_norm.fits'
    psf_norm_path = f'{analysis_folder_path}/{psf_norm_name}'
    
    fits.writeto(psf_norm_path,psf_norm,psf_hdr,overwrite=True)
    
    return psf_norm_name,psf_norm_path


if __name__ == '__main__':

    (analysis_programme,analysis_version,analysis_range,
                    version_continuation,version_control,psf,initial_conditions_mode,
                    pixel_scale,zcal_const,crop_factor) = argparse_values(phase='analysis')


    #--------CODE--------#

    # NEEDED PARAMETERS FOR CREATE A PSF
    fwhm = 2.0336 # FWHM of the PSF
    beta = 1.8781 # Beta parameter of the PSF



    # Obtaining the current working directory
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

    # If PSF is in the folder name a large PSF will be created
    # If no name of PSF is in the folder name a large PSF will be created
    if 'psf' not in folder_name_elements:

        # list to store fits names
        list_fits = []

        # obtain every needed fits from the directory
        for file in sorted(os.listdir(fits_path)):

            # checking for the correct files
            if '.fits' in file and 'psf' not in file:
            
                list_fits.append(file)
                
        # choosing the first fits to obtain the shape
        fits_name = list_fits[0]

        # opening the fits
        hdu = fits.open(f'{fits_path}/{fits_name}')

        hdr = hdu[0].header
        img = (hdu[0].data)[0]

        # obtaining the shape
        x_len = img.shape[0]
        y_len = img.shape[1]

        # galfit needs to change x by y positions
        img_shape = (y_len,x_len)
        img_center = (y_len//2,x_len//2)

        # creating the script for running galfit
        script = f'{folder_new_version_path}/psf_{galaxy}_large.script'
        f = open(script,'w+')

        # the ouput file name will be
        output_file = f'{folder_new_version_path}/psf_{galaxy}_large.fits'

        # calling the function to create the script
        create_psf_script(file=f,
                input_file_name='none',
                output_file_name = output_file,
                psf = 'none',
                cons_file = 'none',
                file_shape = img_shape,
                zp_const = zcal_const,
                pix_scl = pixel_scale,
                img_center = img_center,
                int_mag = 1,
                fwhm = fwhm,
                beta = beta,
                ax_rat = 1,
                pos_ang = 0)

        # closing the script file
        f.close()

        # creating the PSF with the given parameters
        subprocess.run(['galfit', script])

        # opening the result
        hdu = fits.open(f'{output_file}')

        img = hdu[0].data
        hdr = hdu[0].header

        # total flux of the PSF
        total_flux = np.sum(img)

        # normalizing the PSF
        img_norm = img / total_flux

        # saving the result
        fits.writeto(f'{folder_new_version_path}/psf_{galaxy}_large_flux_norm.fits',img_norm,hdr,overwrite=True)


    # If psf is in the folder name a small psf will be created
    else:

        print('\n#------------------IMPORTANT------------------#')
        print('The psf used for the analysis will be the small')
        print(f'Change the folder name to {galaxy}_{version}_PSF\nor to {galaxy}_{version} to create a large PSF for {galaxy}\n')