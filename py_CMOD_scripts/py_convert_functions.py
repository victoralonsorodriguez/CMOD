import numpy as np

from astropy.io import fits

from cmod.io import open_fits
from cmod.utils import round_number


###------------------ELLIPTICITY AND AXIS RATIO------------------###
####################################################################

def ell_to_axrat(ell):
    axrat = 1 - ell
    return round_number(axrat,3)

def axrat_to_ell(axrat):
    ell = 1 - axrat
    return round_number(ell,3)


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


###------------------PIXELS, ARCSEC AND PARSECS------------------###
####################################################################

# Changing from pixeles to kpc 
def px_to_kpc(px,inst):

    arcsec = px * inst 
    kpc = arcsec_to_kpc(arcsec)

    return round_number(kpc,3)

def kpc_to_px(kpc,inst):
    
    arcsec = kpc_to_arcsec(kpc)
    px = arcsec / inst 
    
    return round_number(px,3)

def arcsec_to_kpc(arcsec):
    
    kpc = arcsec * 0.12
    
    return round_number(kpc,3)

def kpc_to_arcsec(kpc):
    
    arcsec = kpc / (0.12)
    
    return round_number(arcsec,3)



###------------------RADIANS AND DEGREES------------------###
#############################################################

# Changing from degrees to radian
def deg_to_rad(deg):
    
    rad = (deg) * (np.pi / 180)

    return round_number(rad,3)

def rad_to_deg(rad):
    
    deg = (rad) * (180 / np.pi)
        
    return round_number(deg,3)


def rad_to_deg_abs(rad):
    
    deg = ((rad) * (180 / np.pi)%360)
        
    return round_number(deg,3)


###------------------REDSHIFT AND WAVELENGHT------------------###
#################################################################

# Transforming the redshift into wavelength
def RtW(z,filter_wl):

    wavelength = (filter_wl)/(1+z)
    
    np.seterr(divide = 'ignore') 
    np.seterr(invalid = 'ignore')
    
    return round_number(wavelength,2)

# Transforming the wavelength into redshift
def WtR(wl,filter_wl):
    
    redshift = ((filter_wl)/wl) - 1 

    np.seterr(divide = 'ignore') 
    np.seterr(invalid = 'ignore')
    
    return round_number(redshift,2)