
import numpy as np
import sys
import os
import subprocess

from astropy.io import fits

from scipy.optimize import curve_fit

from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture
from photutils.isophote import Ellipse

from .fitting import create_psf_script
from .io import open_fits, argparse_values
from .plotting import plot_images, plot_profiles, plot_gen
from .utils import round_number


def isophote_fitting(img_gal_path,
                     gal_center,
                     gal_size,
                     img_mask_path=None,
                     cons=None,
                     output_path='.'):
    
    print('\n\tPerforming the isophote fitting\n')
    
    # Loading the galaxy image
    _,gal_img,gal_img_name = open_fits(img_gal_path)
        
    gal_img_fit = gal_img
    
    # Loading the mask if it is required
    if img_mask_path != None:
        # Loading the mas
        _,mask_img,_ = open_fits(img_mask_path)
        # Converting from 0/1 mask to True/False mask
        mask_img = mask_img == 1
    
        # Masked galaxy image
        img_gal_mask = np.ma.masked_array(gal_img, mask=mask_img)
        gal_img_fit = img_gal_mask
    
    fig_name = f'{gal_img_name}_image_analyze'
    plot_images(gal_img_fit,fig_name,cons=cons,counts=True,log_scale=True,
                output_path=output_path)
    
    
    # While loop beacuse we don't want to fix ellipse initial conditions
    # So the loop is active until the programe converges
    isophote_table = []
    pa_ind = 0
    pa_ini = 2
    pa_fin = 358
    pa_steps = (pa_fin - pa_ini + 1)
    pa_range = np.linspace(pa_ini,pa_fin,pa_steps)
    
    gal_img_fit[np.isnan(gal_img_fit)] = 0

    sma = max(gal_size) * 0.2

    # Image should be in counts
    while len(isophote_table) == 0:

        
        # Defining a elliptical geometry
        print(f'\tAtempting to converge with PA={pa_range[pa_ind]:.2f} deg')
        
        geometry = EllipseGeometry(x0=gal_center[1], y0=gal_center[0],
                                   sma=sma, # semimajor axis in pixels
                                   eps=0.5,
                                   pa=pa_range[pa_ind] * np.pi / 180.0) # position angle in radians
        
        # Creating an aperture object
        aperture = EllipticalAperture((geometry.x0, geometry.y0),
                                      geometry.sma,
                                      geometry.sma * (1 - geometry.eps),
                                      geometry.pa)  
        
        # Fitting the isophotes by using ellipses
        #try:
        fit_step = 0.01        
        ellipse = Ellipse(gal_img_fit, geometry)
        isolist = ellipse.fit_image(step=fit_step,
                                    minit=30,
                                    sclip=3,
                                    nclip=10,
                                    fix_center=True,
                                    fflag=0.5
                                    )

        # We can generate a table with the results
        # omiting the first row to avoid 0 values
        isophote_table = isolist[1:].to_table()
            
        #except Exception as e:
        #    print(f'There is an error in fitting isophotes' 
        #          f'with PA={pa_range[pa_ind]:.2f} deg:'
        #          f'{e}')
            
        
        pa_ind += 1

        if pa_ind == len(pa_range):
            print('Isophote fitting cannot converge')
            break
                
    print('Isophote fitting is finished')
                
    # Export it as a csv
    isophote_table_name = f'{gal_img_name}_isophote.csv'
    isophote_table_path = f'{output_path}/{isophote_table_name}'
    isophote_table.write(f'{isophote_table_path}', format='csv', overwrite=True)
    
    # Creating some figures    
    if 'model' not in gal_img_name:
        plot_list = [(isophote_table_path,'Data')]
        plot_profiles(csv_path_list = plot_list,
                      fig_name = f'{gal_img_name}',cons=cons,
                      output_path=output_path)
        
    else:
        plot_list = [(isophote_table_path,'Model')]
        plot_profiles(csv_path_list = plot_list,
                      fig_name = 'i_model',
                      cons=cons,
                      output_path=output_path)
    
    
    fig_name = f'{gal_img_name}_image_ellipses'
    plot_images(gal_img,fig_name,cons=cons,ellip=True,isolist=isolist,log_scale=True,
                output_path=output_path)
    
    return isophote_table_path



def max_center_value(data,
                     crop_factor=10):


    # Obtaining the shape of the frame
    img_galfit_shape = data.shape
    nrow = data.shape[0]
    ncol = data.shape[1]
    
    if nrow<100 or ncol<100:
        max_pix_value_center = np.nanmax(data)
    else:
        # We are going to crop the frame as we want to obtain the center
        # that will correspond with the maximun intensity value in its center
        # The crop is to avoid external stars
        nrow_range = (nrow//crop_factor)
        ncol_range = (ncol//crop_factor)
        
        # Obtaining the nearest to the center maximun pixel value position to center the model
        frame_center_area = data[nrow//2 - nrow_range:nrow//2 + nrow_range, ncol//2 - ncol_range:ncol//2 + ncol_range]
        
        max_pix_value_center = np.nanmax(frame_center_area)
        
        
    max_pix_value_pos_x = np.where(data == max_pix_value_center)[0][0]
    max_pix_value_pos_y = np.where(data == max_pix_value_center)[1][0]
    
    # In ds9 these values should be +1
    max_pix_value_center_pos = (max_pix_value_pos_x ,max_pix_value_pos_y)
    
    return round_number(max_pix_value_center,3), max_pix_value_center_pos


def sersic_profile(r, n, R_e, I_e):

    b_n = 2 * n - 1 / 3 + 0.009876 / n
    I_r = I_e * np.exp(-b_n * ((r / R_e) ** (1 / n) - 1))
    
    return I_r

def sersic_profile_log(r, n, R_e, a):
    if 0.6<n<6 and 1<R_e<100 and 1<a<6.9:
        b_n = 2 * n - 1 / 3 + 0.009876 / n
        I_r_log = a + (-b_n * ((r / R_e) ** (1 / n) - 1))
        
    elif 0.6<n<6 and 1<R_e<100:
        b_n = 2 * n - 1 / 3 + 0.009876 / n
        I_r_log = a + (-b_n * ((r / R_e) ** (1 / n) - 1))

    else:
        I_r_log = 10e10

    return I_r_log

def exponential_disk(r, I_0, h):
    I_r = I_0 * np.exp(-r / h)
    return I_r

def exponential_disk_linear(x, m, b):
    y = m * x + b
    return y

def two_line_model(x, a1, b1, a2, b2, x_break):
    
    return np.where(x < x_break, a1 * x + b1, a2 * x + b2)

def calculate_chi_squared(x, y_obs, y_err, model, params):

    y_model = model(x, *params)  
    residuals = (y_obs - y_model) / y_err  
    chi_squared = np.sum(residuals**2)  
    return chi_squared


def initial_conditions(img_name,df, x_col, y_col,y_err_col,const,output_path):
    
    print('Computing the initial conditions')
    
    # Loading original data
    x_data = df[x_col].values
    y_data = df[y_col].values
    
    # FIRST STEP: Fitting two straight lines
    # The break point is computed automatically
    print(f'\nFittig two straight lines')
    
    # Converting original data into ln to work with
    # straight lines for a disk like part
    x_data_two = x_data
    y_data_two = np.log(df[y_col].values)
    y_err_two = np.log(df[y_err_col].values)
    
    # Initial estimation for parameters: a1, b1, a2, b2, x_break
    initial_guess_two = [1, 0, -1, 0, max(x_data_two) * (1/3)]
    
    # Fitting the two straight lines model
    popt_two, pcov_two = curve_fit(two_line_model, x_data_two, y_data_two, 
                                   p0=initial_guess_two,sigma=y_err_two)
    
    a1, b1, a2, b2, x_break_two = popt_two
    
    # Computing the ordinary parameter form
    I_0_two = np.exp(b2)
    h_two = -1 / a2
    
    print(f'I_0: {I_0_two:.3f}')
    print(f'h: {h_two:.3f}')
    print(f'X break point: {int(x_break_two)}')

    chi_squared_two = calculate_chi_squared(x_data_two, y_data_two, y_err_two, two_line_model, popt_two)
    print(f"Chi^2: {chi_squared_two:.3f}")
    
    # Obtainting the break position between the two lines
    break_pos = list(x_data_two).index(min(x_data_two[x_data_two >= x_break_two]))
    print(f'X break point pos: {break_pos}')

    # Representing the lines
    y_fit_line_1 = a1 * x_data_two[x_data_two < x_break_two] + b1
    y_fit_line_2 = a2 * x_data_two[x_data_two >= x_break_two] + b2
    
    y_fit_line_1 = list(y_fit_line_1) + [np.nan] * (len(x_data_two) - len(y_fit_line_1))
    y_fit_line_2 = [np.nan] * (len(x_data_two) - len(y_fit_line_2)) + list(y_fit_line_2)
    
    x_data_fig = [x_data,x_data,x_data]
    y_data_fig = [y_data_two,y_fit_line_1,y_fit_line_2]
    
    # Plotting the lines
    fig_name = f'{img_name}_IC_00'
    label_list = ['Intensity Profile','Line 1', 'Line 2']
    
    plot_gen(x_data=x_data_fig,
             y_data=y_data_fig,
             label_list=label_list,
             fig_name=fig_name,
             fig_save_path=output_path,
             line_color = ['gold','dodgerblue','black'],
             x_axis_label = ['$X\\, dimension\\, [\\mathrm{{pix}}]$'],
             x_ticks_dec = 0,
             y_axis_label = ['$\\ln(Intensity/\\mathrm{{counts}})$'],
             y_ticks_dec = 2,
             width_style = [0.75,2,2],
             zorder = [1,2,2],
             style_plot = ['sct','line','line'],
             plot_show= False)
    
    
    # SECOND STEP: disk fitting
    print(f'\nDisk-like fitting')
    
    # The disk-like fitting will be carry out with the
    # external part from the break position
    # Transforming data to ln
    x_data_disk = df[x_col].values[break_pos:]
    y_data_disk = np.log(df[y_col].values[break_pos:])
    y_err_disk = np.log(df[y_err_col].values[break_pos:])

    initial_guess_disk = [-1, np.mean(x_data_disk)]
    
    # Fitting the model to an exponential converted to
    # a straight line
    popt_disk, pcov_disk = curve_fit(exponential_disk_linear, x_data_disk, y_data_disk, 
                                     p0=initial_guess_disk,sigma=y_err_disk,nan_policy='omit')
    m_disk, b_disk = popt_disk
    
    # Changing from the straight line to
    # the exponential parameters
    I_0_disk = np.exp(b_disk)
    h_disk = -1/m_disk

    print(f'I_0: {I_0_disk:.3f}')
    print(f'h: {h_disk:.3f}')
        
    chi_squared_disk = calculate_chi_squared(x_data_disk, y_data_disk, y_err_disk, exponential_disk_linear, popt_disk)
    print(f"Chi^2 disk: {chi_squared_disk:.3f}")
    
    # Creating the exponential 
    y_fit_disk = exponential_disk_linear(x_data_disk, m_disk, b_disk)
    y_fit_disk = exponential_disk(x_data_disk,I_0_disk,h_disk)
    
    
    # THIRD STEP: bulge fitting
    print(f'\nBulge-like fitting')

    # For this fit an exponential is applied
    # affecting only to the intensity
    break_pos_bul = int(break_pos)
    x_data_bul = df[x_col].values[:break_pos_bul]
    y_data_bul = np.log(df[y_col].values[:break_pos_bul])
    y_err_bul = np.log(df[y_err_col].values[:break_pos_bul])
    
    if '-hr' in sys.argv:
        
        x_data_bul = [val for val in x_data_bul if list(x_data_bul).index(val)%2 == 0]
        y_data_bul = [val for val in y_data_bul if list(y_data_bul).index(val)%2 == 0]
        y_err_bul = [val for val in y_err_bul if list(y_err_bul).index(val)%2 == 0]

    initial_guess_bul = [3,40,3]
    
    # Fitting the model
    popt_bul, pcov_bul = curve_fit(sersic_profile_log, x_data_bul, y_data_bul, 
                                   p0=initial_guess_bul,sigma=y_err_bul)
    n_bul, r_e_bul, a = popt_bul
    
    # Computing the ordinary parameter form
    I_e_bul = np.exp(a)
    
    print(f'n: {n_bul:.3f}')
    print(f'R_e: {r_e_bul:.3f}')
    print(f'I_e: {I_e_bul:.3f}')
    chi_squared_bul = calculate_chi_squared(x_data_bul, y_data_bul, y_err_bul, sersic_profile_log, popt_bul)
    print(f'Chi^2: {chi_squared_bul:.3f}')
    
    # Crear los ajustes
    y_fit_bul = sersic_profile(x_data_bul, n_bul, r_e_bul, I_e_bul)
    
    # Ploting the fitting
    
    y_fit_bul = np.array(list(y_fit_bul) + [np.nan] * (len(x_data_two) - len(y_fit_bul)))
    y_fit_disk = np.array([np.nan] * (len(x_data_two) - len(y_fit_disk)) + list(y_fit_disk))
    
    x_data_fig = [x_data,x_data,x_data]
    y_data_fig = [y_data,y_fit_bul,y_fit_disk]
    
    label_list = ['Intensity Profile','Sersic Fitting - Bulge-like', 'Exponential Fitting - Disk-like']
    fig_name = f'{img_name}_IC_01'    
    
    plot_gen(x_data=x_data_fig,
             y_data=y_data_fig,
             label_list=label_list,
             fig_name=fig_name,
             fig_save_path=output_path,
             line_color = ['gold','dodgerblue','black'],
             x_axis_label =['$X\\, dimension\\, [\\mathrm{{pix}}]$'],
             x_ticks_dec = 0,
             y_axis_label = ['$Intensisty\\, [\mathrm{{counts}}]$'],
             y_ticks_dec = 2,
             width_style = [0.75,2,2],
             zorder = [1,2,2],
             style_plot = ['sct','line','line'],
             plot_show=False)
    
    return (break_pos, break_pos_bul, round_number(I_0_disk,3), round_number(h_disk,3), 
            round_number(n_bul,3), round_number(r_e_bul,3), round_number(I_e_bul,3))
    

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