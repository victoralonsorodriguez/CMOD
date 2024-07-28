from astropy.io import fits
import numpy as np

import os
import subprocess

from py_psf_script import create_psf_script

# This script creates a large PSF with the
# size of the galaxy image

# NEEDED DATA
cal = 25 # zeropoint calibration constant
inst = 0.2 # pixel scale arcsec / pix
fwhm = 2.0336 # FWHM of the PSF
beta = 1.8781 # Beta parameter of the PSF


# obtaining the current working directory
cwd = os.getcwd()
galaxy = (cwd.split('/')[-1]).split('_')[0]


# list to store fits names
list_fits = []

# obtain every needed fits from the directory
for file in sorted(os.listdir(cwd)):

	# checking for the correct files
	if '.fits' in file and 'psf' not in file:
	
	    list_fits.append(file)
	    
# choosing the first fits to obtain the shape
fits_name = list_fits[0]

# opening the fits
hdu = fits.open(f'{cwd}/{fits_name}')

hdr = hdu[0].header
img = (hdu[0].data)[0]

# obtaining the shape
x_len = img.shape[0]
y_len = img.shape[1]

# galfit needs to change x by y positions
img_shape = (y_len,x_len)
img_center = (y_len//2,x_len//2)

# creating the script for running galfit
script = f'psf_{galaxy}.script'
f = open(script,'w+')

# the ouput file name will be
output_file = f'psf_{galaxy}_large.fits'

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
hdu = fits.open(f'{cwd}/{output_file}')

img = hdu[0].data
hdr = hdu[0].header

# total flux of the PSF
total_flux = np.sum(img)

# normalizing the PSF
img_norm = img / total_flux

# saving the result
fits.writeto(f'{cwd}/psf_{galaxy}_large_flux_norm.fits',img_norm,hdr,overwrite=True)


