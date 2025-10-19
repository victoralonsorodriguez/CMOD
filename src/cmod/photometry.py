import numpy as np
from astropy.io import fits

from .io import open_fits
from .utils import round_number

###------------------MAGNITUDES AND COUNTS------------------###
###############################################################

def fits_mag_to_counts(data = None,
                       fits_input_path = None,
                       fits_output_path = '.',
                       const = (0.1,20),
                       remove_padding = False,
                       padding_value=np.nan):
    
    if data == None and fits_input_path != None:
        hdr_mag,data_mag,fits_mag_name = open_fits(fits_input_path)
    elif data != None:
        data_mag = data
    else:
        print('Data is requiered to convert magnitudes into counts')

    # Flux from magnitudes to counts
    fits_flux_img = 10**((data_mag - 2.5*np.log10(const[0]**2)-const[1])/-2.5)
    
    # Deleting the padding 
    if remove_padding == True:
        padding = fits_flux_img == fits_flux_img[0,0]
        fits_flux_img[padding] = padding_value
    
    # Expoting the new fits
    if 'mag' in fits_mag_name:
        fits_flux_name = fits_mag_name.replace('mag', 'counts')+'.fits'
    else:
        fits_flux_name = f'{fits_mag_name}_counts.fits'

    fits_flux_path = f'{fits_output_path}/{fits_flux_name}'
    fits.writeto(f'{fits_flux_path}', fits_flux_img, header=hdr_mag, overwrite=True)
    
    return fits_flux_path



def fits_counts_to_mag(data = None,
                       fits_input_path = None,
                       fits_output_path = '.',
                       const = (0.1,20),
                       remove_padding = False,
                       padding_value=np.nan):
    
    if data == None and fits_input_path != None:
        hdr_flux,data_flux,fits_flux_name = open_fits(fits_input_path)
    elif data != None:  
        data_flux = data
    else:
        print('Data is requiered to convert magnitudes into counts')

    # Flux from counts to mag/arcsec^2
    data_flux[data_flux<=0] = 1
    fits_mag_img = -2.5*np.log10(data_flux) + 2.5 * np.log10(const[0]**2) + const[1]
    
    # Deleting the padding 
    if remove_padding == True:
        padding = fits_mag_img == fits_mag_img[0,0]
        fits_mag_img[padding] = padding_value
    
    # Expoting the new fits in mag
    if 'counts' in fits_flux_name:
        fits_flux_name = fits_flux_name.replace('counts', 'mag')+'.fits'
    else:
        fits_flux_name = f'{fits_flux_name}_mag.fits'
        
    fits_mag_path = f'{fits_output_path}/{fits_flux_name}'
    fits.writeto(f'{fits_mag_path}', fits_mag_img, header=hdr_flux, overwrite=True)
    
    return fits_mag_path  


def values_counts_to_mag(val_counts,const = (0.1,20)):
    
    if isinstance(val_counts, (list,np.ndarray)):
        
        val_counts[val_counts<=0]=1
        val_mag = -2.5*np.log10(val_counts) + const[1]
        
        return round_number(val_mag,2)
        
    elif isinstance(val_counts, (int, float, np.float32, np.float64)) and (val_counts<=0 or np.isnan(val_counts)):

        val_counts = 1
        val_mag = -2.5*np.log10(val_counts) + const[1]
        
        return round_number(val_mag,2)
    
    elif isinstance(val_counts, (int, float, np.float32, np.float64)) and val_counts>0:

        val_mag = -2.5*np.log10(val_counts) + const[1]
        
        return round_number(val_mag,2)
    
    else:
        print(val_counts)
        print(np.isnan(val_counts))
        print(type(val_counts))
        print('No valid value')


def values_mag_to_counts(val_mag,const):
    
    val_count = 10**((val_mag - 2.5*np.log10(const[0]**2)-const[1])/-2.5)
    
    return round_number(val_count,3)