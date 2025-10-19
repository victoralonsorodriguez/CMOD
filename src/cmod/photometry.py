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


def filter_wl_dict():
    
    # The names of the filter systems in alphabetic order
    filter_system_names = ['Euclid','JWST','Hubble','Johnson_UVRIJHK']
    filter_system_corrected_names = ['Euclid','JWST','HST','SDSS']
    
    # Dictionary with the corresponding wavelength of each filter
    # http://svo2.cab.inta-csic.es/theory/fps/
    filter_wavelenght_dict = {'EucHab':14971.70, 'EucJab':11522.58 , 'EucVISab':4970.77,
                'EucYab':9381.52,
                'F070W':7039.12, 'F090W':9021.53, 'F115W':11542.61, 
                'F140M':14053.23, 'F150W':15007.44, 'F162M':16272.47, 
                'F182M':18451.67, 'F200W':19886.48, 'F210M':20954.51, 
                'F277W':27617.40, 'F356W':35683.62, 'F444W':44043.15,
                'Hve':16396.38, 'Ive':8657.44, 'Jve':12317.30, 'Kve':21735.85,
                'Rve':6819.05, 'Uve':3511.89, 'Vve':5501.40,
                'HSTF300W':2984.47, 'HSTF435W':4323.09, 'HSTF450W':4556.22,
                'HSTF475W':4775.73, 'HSTF555W':5355.74, 'HSTF569W':5644.08, 
                'HSTF606W':5887.08, 'HSTF791W':7875.56, 
                'HSTF814W':8039.03}
    
    # Dictionary with the correct name to plot for each filter
    filter_names_dict = {'EucHab':'HE', 'EucJab':'JE' , 'EucVISab':'IE',
                'EucYab':'YE',
                'F070W':'F070W', 'F090W':'F090W', 'F115W':'F115W', 
                'F140M':'F140M', 'F150W':'F150W', 'F162M':'F162M', 
                'F182M':'F182M', 'F200W':'F200W', 'F210M':'F210M', 
                'F277W':'F277W', 'F356W':'F356W', 'F444W':'F444W',
                'Hve':'H', 'Ive':'I', 'Jve':'J', 'Kve':'K',
                'Rve':'R', 'Uve':'U', 'Vve':'V',
                'HSTF300W':'F300W', 'HSTF435W':'F435W', 'HSTF450W':'F450W',
                'HSTF475W':'F475W', 'HSTF555W':'F555W', 'HSTF569W':'F569W', 
                'HSTF606W':'F606W', 'HSTF791W':'F791W', 
                'HSTF814W':'F814W'}
                
                
    return filter_system_corrected_names,filter_wavelenght_dict,filter_names_dict