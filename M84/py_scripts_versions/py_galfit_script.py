
# This file creates the required template for Galfit
def create_script(file,input_file_name,output_file_name,psf,cons_file,file_shape,zp_const,pix_scl,img_center,int_mag,eff_rad,ser_ind,ax_rat,pos_ang,noise_file = "none",pixel_mask= "none"):
	
	file.write("# IMAGE PARAMETERS\n")
	file.write(f" A) {input_file_name} # Input Data image (FITS file))\n")
	file.write(f" B) {output_file_name} # Name for the output image)\n")
	file.write(f" C) {noise_file} # Noise image name (made from data if blank or 'none')\n")
	file.write(f" D) {psf} # Input PSF image and (optional) diffusion kernel\n")
	file.write(f" E) 1 # PSF oversampling factor relative to data\n")
	file.write(f" F) {pixel_mask} # Pixel mask (ASCII file or FITS file with non-0 values)\n")
	file.write(f" G) {cons_file} # Parameter constraint file (ASCII)\n")
	file.write(f" H) 1		{file_shape[1]} 1		{file_shape[0]} # Image region to fit (xmin xmax ymin ymax)\n")
	file.write(f" I) {file_shape[1]}		{file_shape[0]} # Size of convolution box (x y)\n")
	file.write(f" J) {zp_const} # Magnitude photometric zeropoint\n")
	file.write(f" K) {pix_scl} {pix_scl} # Plate/pixel scale (dx dy)\n")
	file.write(f" O) regular # Display type (regular, curses, both)\n")
	file.write(f" P) 0 # Create ouput only? (1=yes; 0=optimize)\n")
	file.write(f" S) 0 # Modify/create objects interactively?\n")
	
	file.write(f"\n# Component number: 1\n")
	file.write(f" 0) sersic                 #  Component type\n")
	file.write(f" 1) {img_center[1]} {img_center[0]} 1 1 #  Position centre x, y\n")
	file.write(f" 3) {int_mag}     1          #  Integrated magnitude\n")
	file.write(f" 4) {eff_rad}     1          #  R_e (effective radius)   [pix]\n")
	file.write(f" 5) {ser_ind}     1          #  Sersic index n (de Vaucouleurs n=4)\n")
	file.write(f" 9) {ax_rat}      1          #  Axis ratio (b/a)\n")
	file.write(f"10) {pos_ang}     1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
	file.write(f" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
	
	
