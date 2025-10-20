# En src/cmod/resampling.py
import numpy as np
from scipy.ndimage import zoom
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import cv2
from skimage.transform import resize

from cmod.cosmology import zlocal_correction, cosmological_scale
from cmod.io import open_fits

def resampling_frame_pixels(galaxy,
                            galaxy_shape,
                            original_sensor_pixscale,
                            zsim,
                            simulated_sensor_pxscale):
    
    # Local galaxy reshift and scale
    galaxy_zlocal = zlocal_correction(galaxy)
    galaxy_scale = cosmological_scale(galaxy_zlocal)
    
    # Lcoal galaxy pixel scale and total field of view
    galaxy_pixscale = galaxy_scale * original_sensor_pixscale
    galaxy_fov_X = galaxy_pixscale * galaxy_shape[1]
    galaxy_fov_Y = galaxy_pixscale * galaxy_shape[0]
    
    # For a simulated redshift: pixel scale and field of view
    if zsim == 0.0:
        zsim = galaxy_zlocal
    zsim_scale = cosmological_scale(zsim)
    
    zsim_fov_X = galaxy_fov_X / zsim_scale
    zsim_fov_Y = galaxy_fov_Y / zsim_scale
    
    # Computing pixel side of resampled image
    pixel_size_X = int(round(zsim_fov_X / simulated_sensor_pxscale,0))
    pixel_size_Y = int(round(zsim_fov_Y / simulated_sensor_pxscale,0))
    
    resample_frame_size = (pixel_size_Y,pixel_size_X)
    
    
    return resample_frame_size



def resampling_frame(frame_path,resample_frame_size):
    
    # Loading the fits file
    hdr,img,fits_name = open_fits(frame_path)
    
    print(f'Resampling file: {fits_name}')
    
    # Managing nan values
    nan_mask = np.isnan(img)
    image_no_nan = np.nan_to_num(img, nan=0.0)
    
    # Apliying the zoom function
    resampled_frame_no_nan = resize(image_no_nan,
                                    resample_frame_size,
                                    order=3, 
                                    anti_aliasing=False, 
                                    mode='constant', cval=0.0)
    
    resampled_nan_mask = resize(nan_mask.astype(float),
                                resample_frame_size,
                                order=0, 
                                anti_aliasing=False,
                                mode='constant', cval=0.0)
    
    resampled_frame = resampled_frame_no_nan
    resampled_frame[resampled_nan_mask > 0.5] = np.nan
    
    # Renaming the frame
    resampled_frame_name = f'{fits_name}_r{int(resample_frame_size[1])}p.fits'
    resampled_frame_path = f'{frame_path.split(fits_name)[0]}/{resampled_frame_name}'
    
    fits.writeto(resampled_frame_path, resampled_frame, header=hdr, overwrite=True)
    return resampled_frame_path


if __name__ == '__main__':
    
    galaxy_name = 'M84'
    pixel_scale = 0.2
    sensor_pix = 300
    z = 0.0
    sim_pix_scale = 0.11
    
    pixels,zoom_factor = resampling_frame_pixels(galaxy_name,pixel_scale,sensor_pix,
                                         z,sim_pix_scale)
    
    fits_path = '/Users/victor/TFG/galaxias/CMOD/data/M84_analysis_V10_Galfit/m84_VBIN018_SL_zSimJ_EucHab/m84_VBIN018_SL_zSimJ_EucHab_z0.00_counts.fits' 
    
    resampled_frame_path = resampling_frame(fits_path,zoom_factor,pixels)
