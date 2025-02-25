
import numpy as np
from py_round_number import round_number

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