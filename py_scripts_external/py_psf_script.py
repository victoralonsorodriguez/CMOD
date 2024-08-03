

# This file creates the required template for Galfit to create a PSF based on a Moffat function
def create_psf_script(file,input_file_name,output_file_name,psf,cons_file,file_shape,zp_const,pix_scl,img_center,int_mag,fwhm,beta,ax_rat,pos_ang,noise_file = 'none',pixel_mask= 'none'):
	
	file.write("# IMAGE PARAMETERS and GALFIT CONTROL PARAMETERS\n")
	file.write(f" A) {input_file_name} # Input Data image (FITS file))\n")
	file.write(f" B) {output_file_name} # Name for the output image)\n")
	file.write(f" C) {noise_file} # Noise image name (made from data if blank or 'none')\n")
	file.write(f" D) {psf} # Input PSF image and (optional) diffusion kernel\n")
	file.write(f" E) 1 # PSF oversampling factor relative to data\n")
	file.write(f" F) {pixel_mask} # Pixel mask (ASCII file or FITS file with non-0 values)\n")
	file.write(f" G) {cons_file} # Parameter constraint file (ASCII)\n")
	file.write(f" H) 1		{file_shape[0]} 1		{file_shape[1]} # Image region to fit (xmin xmax ymin ymax)\n")
	file.write(f" I) {file_shape[0]}		{file_shape[1]} # Size of convolution box (x y)\n")
	file.write(f" J) {zp_const} # Magnitude photometric zeropoint\n")
	file.write(f" K) {pix_scl} {pix_scl} # Plate/pixel scale (dx dy)\n")
	file.write(f" O) regular # Display type (regular, curses, both)\n")
	file.write(f" P) 1 # Options: 0=normal run; 1,2=make model/imgblock & quit\n")
	
	file.write(f"\n# Moffat function\n")
	file.write(f" 0) moffat     # object type\n")
	file.write(f" 1) {img_center[0]} {img_center[1]} 1 1 #  Position centre x, y [pix]\n")
	file.write(f" 3) {int_mag}     1          #  total magnitude\n")
	file.write(f" 4) {fwhm}     1          #  FWHM   [pix]\n")
	file.write(f" 5) {beta}     1          #  powerlaw (beta)\n")
	file.write(f" 9) {ax_rat}      1          #  axis ratio (b/a)   \n")
	file.write(f"10) {pos_ang}     1          #  position angle (PA)  [Degrees: Up=0, Left=90] \n")
	file.write(f" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
	
	
