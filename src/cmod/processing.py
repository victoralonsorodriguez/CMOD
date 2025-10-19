
import numpy as np

from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture
from photutils.isophote import Ellipse

from .io import open_fits
from .plotting import plot_images, plot_profiles
from .utils import round_number


def isophote_fitting(img_gal_path,
                     gal_center,
                     img_mask_path=None,
                     cons=None,
                     output_path='.'):
    
    print('Performing the isophote fitting\n')
    
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
    pa_fin = 178
    pa_steps = (pa_fin - pa_ini + 1)
    pa_range = np.linspace(pa_ini,pa_fin,pa_steps)
    
    
    gal_img_fit[np.isnan(gal_img_fit)] = 0

    iso_fitted = False

    # Image should be in counts
    #while iso_fitted == False:
    while len(isophote_table) == 0:

        
        # Defining a elliptical geometry
        print(f'Atempting to converge with PA={pa_range[pa_ind]:.2f} deg')
        
        geometry = EllipseGeometry(x0=gal_center[1], y0=gal_center[0],
                                   sma=50, # semimajor axis in pixels
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