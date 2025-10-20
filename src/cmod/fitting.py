


# Script to store Galfit constraints values
def constraints_values():

    n_range = (0.6,6)
    re_range = (1,600)
    mag_range = (-1.0,1.0)
    xy_range = (-2.0,2.0)
    q_range = (0,1)
    pa_range = (-90.00,90.00)

    return n_range,re_range,mag_range,xy_range,q_range,pa_range


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
 
 
# This file creates the required constraints template for Galfit
def create_constraints(file,n_range,re_range,mag_range,xy_range,q_range,pa_range):
	
	file.write('# Component   parameter   constraint  (Optional and written with #)Commentary\n')
	file.write(f'# operation values \n')
	file.write(f'1 n {n_range[0]} to {n_range[1]} #polis email correction \n')
	file.write(f'1 re {re_range[0]} to {re_range[1]}  \n')
	file.write(f'1 mag {mag_range[0]} {mag_range[1]}  \n')
	file.write(f'1 x {xy_range[0]} {xy_range[1]}  \n')
	file.write(f'1 y {xy_range[0]} {xy_range[1]}  \n')
	file.write(f'1 q {q_range[0]} to {q_range[1]} \n')
	file.write(f'1 pa {pa_range[0]} to {pa_range[1]}')




# This file creates the required template for Galfit
def create_conf_imfit(file,funct,galaxy,img_center_x,img_center_y,
                     pa_ser=[None,None,None],ell_ser=[None,None,None],n=[None,None,None],I_e=[None,None,None],r_e=[None,None,None], # Sersic
                     pa_disk=[None,None,None],ell_disk=[None,None,None],I_0=[None,None,None],h=[None,None,None],                      # Exponential
                     pa_bar=[None,None,None],ell_bar=[None,None,None],n_bar=[None,None,None],a_bar=[None,None,None],c0=[None,None,None]): # Ferrer
	
    if len(funct) != 0:
        
        file.write(f'# Configuration file for the galaxy {galaxy}\n')
        file.write(f'X0    {img_center_x[0]}    {img_center_x[1]},{img_center_x[2]}\n')
        file.write(f'Y0    {img_center_y[0]}    {img_center_y[1]},{img_center_y[2]}\n')

        for f in funct:
            if f == 'Sersic':
                file.write(f'FUNCTION {f}\n')
                file.write(f'PA    {pa_ser[0]}    {pa_ser[1]},{pa_ser[2]}\n')
                file.write(f'ell    {ell_ser[0]}    {ell_ser[1]},{ell_ser[2]}\n')
                file.write(f'n    {n[0]}    {n[1]},{n[2]}\n')
                file.write(f'I_e    {I_e[0]}    {I_e[1]},{I_e[2]}\n')
                file.write(f'r_e    {r_e[0]}    {r_e[1]},{r_e[2]}\n')

            elif f == 'Exponential':
                file.write(f'FUNCTION {f}\n')
                file.write(f'PA    {pa_disk[0]}    {pa_disk[1]},{pa_disk[2]}\n')
                file.write(f'ell    {ell_disk[0]}    {ell_disk[1]},{ell_disk[2]}\n')
                file.write(f'I_0    {I_0[0]}    {I_0[1]},{I_0[2]}\n')
                file.write(f'h    {h[0]}    {h[1]},{h[2]}\n')
                
            elif f == 'FerrersBar2D':
                file.write(f'FUNCTION {f}\n')
                file.write(f'PA    {pa_bar[0]}    {pa_bar[1]},{pa_bar[2]}\n')
                file.write(f'ell    {ell_bar[0]}    {ell_bar[1]},{ell_bar[2]}\n')
                file.write(f'I_0    {I_0[0]}    {I_0[1]},{I_0[2]}\n')
                file.write(f'n    {n_bar[0]}    {n_bar[1]},{n_bar[2]}\n')
                file.write(f'a_bar    {a_bar[0]}    {a_bar[1]},{a_bar[2]}\n')
                file.write(f'c0    {c0[0]}    {c0[1]},{c0[2]}\n')

    else:
        print('No functions for the fitting were selected')
        return None
	
	


# This file creates the required constraints template for Galfit
def create_initiaL_params(file,galaxy,version,eff_rad,ser_ind,ax_rat,pos_ang):
	
	file.write(f'# {galaxy} initial parameters introduced manually for {version}\n')
	file.write(f'Effective radius: {eff_rad} [pix]\n')
	file.write(f'Sérsic Index: {ser_ind}  \n')
	file.write(f'Axis ratio: {ax_rat}  \n')
	file.write(f'Position angle: {pos_ang} [degrees, up=0, left=90]  \n')





# Script to store Galfit initial parameters
def initial_params(galaxy):

    if galaxy == 'ESO498G05':
    
        eff_rad = 73
        ser_ind = 4.0
        ax_rat = 0.85
        pos_ang = -23   

    elif galaxy == 'IC719':
    
        eff_rad = 43
        ser_ind = 1.7
        ax_rat = 0.26
        pos_ang = 52

    elif galaxy == 'IC2051':
    
        eff_rad = 63
        ser_ind = 5.0
        ax_rat = 0.65
        pos_ang = 65

    elif galaxy == 'M84':
    
        eff_rad = 60
        ser_ind = 3.0
        ax_rat = 0.81
        pos_ang = -52

    elif galaxy == 'NGC0289':
    
        eff_rad = 76
        ser_ind = 3.0
        ax_rat = 0.60
        pos_ang = -60

    elif galaxy == 'NGC307':
    
        eff_rad = 30
        ser_ind = 2.7
        ax_rat = 0.51
        pos_ang = 80

    elif galaxy == 'NGC788':
    
        eff_rad = 55
        ser_ind = 2.0
        ax_rat = 0.52
        pos_ang = -75
 
    elif galaxy == 'NGC1309':
    
        eff_rad = 67
        ser_ind = 3.5
        ax_rat = 0.89
        pos_ang = 3

    elif galaxy == 'NGC1440':
    
        eff_rad = 47
        ser_ind = 4.0
        ax_rat = 0.8
        pos_ang = 45

    elif galaxy == 'NGC1553':
    
        eff_rad = 75
        ser_ind = 4.5
        ax_rat = 0.66
        pos_ang = -26

    elif galaxy == 'NGC3393':
    
        eff_rad = 75
        ser_ind = 4.2
        ax_rat = 0.66
        pos_ang = -26

    elif galaxy == 'NGC3783':
    
        eff_rad = 67
        ser_ind = 4.0
        ax_rat = 0.8
        pos_ang = -25

    elif galaxy == 'NGC4418':
    
        eff_rad = 50
        ser_ind = 3.8
        ax_rat = 0.60
        pos_ang = 58

    elif galaxy == 'NGC5806':
    
        eff_rad = 66
        ser_ind = 5.0
        ax_rat = 0.56
        pos_ang = -10

    elif galaxy == 'NGC6958':
    
        eff_rad = 66
        ser_ind = 3.9
        ax_rat = 0.56
        pos_ang = -10

    return eff_rad,ser_ind,ax_rat,pos_ang



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