import numpy as np
from cmod.utils import round_number


###------------------ELLIPTICITY AND AXIS RATIO------------------###
####################################################################

def ell_to_axrat(ell):
    axrat = 1 - ell
    return round_number(axrat,3)

def axrat_to_ell(axrat):
    ell = 1 - axrat
    return round_number(ell,3)


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