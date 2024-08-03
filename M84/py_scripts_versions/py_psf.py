
# This script creates a PSF of the size of the galaxy image

from astropy.io import fits
import numpy as np

import os
import subprocess

import pdb

from py_psf_script import create_psf_script

#--------CODE--------#

# NEEDED PARAMETERS FOR CREATE A PSF
cal = 25 # zeropoint calibration constant
inst = 0.2 # pixel scale arcsec / pix
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
			zp_const = cal,
			pix_scl = inst,
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