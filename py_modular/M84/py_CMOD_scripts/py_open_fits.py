
'''
This file contains a function to open fits
'''

# Importing Astropy to open fits files
from astropy.io import fits


def open_fits(fits_path,
              dim = 0):
    
    '''
    This function opens a fits file using astropy.fits.open
    
    Arguments:
        fits_path: str
            path to the original galaxy fits
            
    Returns:
        hdr: astropy.header
            header of the fits
            
        img: numpy.array (2D if image or 3D if datacube)
            data extracted from the fits
            
        galaxy_fits_name: str
            name of galaxy fits without the extension name (.fits)
    '''
    
    # Obtaining the name without the extension
    fits_name = fits_path.split('/')[-1].split('.fits')[0]
    
    # Open the fits with Astropy and sxtracting the header and 
    # the data from it
    hdu = fits.open(f'{fits_path}')
    hdr = hdu[dim].header
    img = hdu[dim].data
    
    # Returning the header, the image and the name
    return (hdr,img,fits_name)