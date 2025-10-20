import numpy as np
from astropy.cosmology import FlatLambdaCDM

from .utils import round_number




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


# This scripts avoids the 0 kpc galaxy radius at z=0
def zlocal_correction(galaxy):

    # Dictionary to store the z=0 kpc/'' factor correction
    zlocal_dict = {'ESO498G05':0.192,
                    'IC719':0.006114, # z = 0.006114
                    'IC2051':0.125, 
                    'M84':0.003392, # z = 0.003392
                    'NGC0289':0.096,
                    'NGC307':0.250,
                    'NGC788':0.266,
                    'NGC1309':0.140,
                    'NGC1440':0.105,
                    'NGC1553':0.075,
                    'NGC3393':0.285,
                    'NGC3783':0.226,
                    'NGC4418':0.174,
                    'NGC5806':0.110,
                    'NGC6958':0.176
                    }
                    
    return zlocal_dict[galaxy]


def cosmological_scale(z):
    
    cosmo_model = FlatLambdaCDM(H0=69.6, Om0=0.286) # https://www.astro.ucla.edu/~wright/CosmoCalc.html
    cosmo_scale = 1./cosmo_model.arcsec_per_kpc_proper(z)
    
    return cosmo_scale.value