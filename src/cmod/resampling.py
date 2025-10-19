# En src/cmod/resampling.py
import numpy as np
from scipy.ndimage import zoom
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

import cv2
from skimage.transform import resize

from cmod.cosmology import kpc_correction, cosmological_scale
from cmod.io import open_fits

def resampling_frame_pixels(galaxy,
                            original_sensor_pixscale,
                            original_sensor_pix_num,
                            zsim,
                            simulated_sensor_pxscale):
    
    # Local galaxy reshift and scale
    galaxy_zlocal = kpc_correction(galaxy)
    galaxy_scale = cosmological_scale(galaxy_zlocal)
    
    # Lcoal galaxy pixel scale and total field of view
    galaxy_pixscale = galaxy_scale * original_sensor_pixscale
    galaxy_fov = galaxy_pixscale * original_sensor_pix_num

    # For a simulated redshift: pixel scale and field of view
    if zsim == 0.0:
        zsim = galaxy_zlocal
    zsim_scale = cosmological_scale(zsim)
    zsim_fov = galaxy_fov / zsim_scale
    
    # Computing pixel side of resampled image
    resample_frame_pixels = int(round(zsim_fov / simulated_sensor_pxscale,0))
    
    # Computing zoom factor
    zoom_factor = resample_frame_pixels / original_sensor_pix_num 
    
    return resample_frame_pixels



def resampling_frame(frame_path,resample_frame_pixels):
    
    # Loading the fits file
    hdr,img,fits_name = open_fits(frame_path)
    
    print(f'Resampling file: {fits_name}')
    
    # Managing nan values
    nan_mask = np.isnan(img)
    image_no_nan = np.nan_to_num(img, nan=0.0)
    
    # Apliying the zoom function
    resampled_frame_no_nan = resize(image_no_nan,
                                    (resample_frame_pixels,resample_frame_pixels),
                                    order=3, 
                                    anti_aliasing=False, 
                                    mode='constant', cval=0.0)
    
    resampled_nan_mask = resize(nan_mask.astype(float), # Convertir a float para resize
                             (resample_frame_pixels,resample_frame_pixels),
                             order=0, 
                             anti_aliasing=False,
                             mode='constant', cval=0.0)
    
    resampled_frame = resampled_frame_no_nan
    resampled_frame[resampled_nan_mask > 0.5] = np.nan
    
    # Renaming the frame
    resampled_frame_name = f'{fits_name}_r{int(resample_frame_pixels)}p.fits'
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
