import numpy as np
from cmod.utils import round_number


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