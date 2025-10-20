
import pandas as pd
import numpy as np
import re

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

from .cosmology import zlocal_correction, cosmological_scale
from .io import open_fits
from .utils import round_number



def galfit_init_dataframe():
    
    # creating the Pandas dataframe to store the information
    df = pd.DataFrame(columns = ['z',
                                 'X_cent', 
                                 'X_cent_err', 
                                 'Y_cent',
                                 'Y_cent_err', 
                                 'Mag','Mag_err',
                                 'Eff_rad',
                                 'Eff_rad_err',
                                 'Eff_rad_arcsec',
                                 'Eff_rad_arcsec_err',
                                 'zscale',
                                 'Eff_rad_kpc',
                                 'Eff_rad_kpc_err',
                                 'Ind','Ind_err',
                                 'Ax_rat',
                                 'Ax_rat_err',
                                 'Pos_ang',
                                 'Pos_ang_err',
                                 'Chi2',
                                 'Chi2nu'])

    return df

"""
def galfit_create_dataframe(df,galaxy_name,output_file_path,pxsc_zcal_const):
    
    print('Creating a dataframe to store results')
    
    hdr,img,frame_output_name = open_fits(output_file_path,dim=2)
    
    # find the redshif in the ouput file name
    z = re.findall('z\d.\d\d',output_file_path)[0][1:]
    
    # output values from galfit and their errors
    x_center = hdr['1_XC'].split(' ')[0]
    x_center_err = hdr['1_XC'].split(' ')[2]
    
    y_center = hdr['1_YC'].split(' ')[0]
    y_center_err = hdr['1_YC'].split(' ')[2]
    
    mag = hdr['1_MAG'].split(' ')[0]
    mag_err = hdr['1_MAG'].split(' ')[2]

    if '*' in mag:
        mag = np.nan
        mag_err = np.nan
    
    eff_rad = hdr['1_RE'].split(' ')[0]
    eff_rad_err = hdr['1_RE'].split(' ')[2]

    if '*' in eff_rad:
        eff_rad = np.nan
        eff_rad_err = np.nan
    
    # Effective radius in arcsec
    eff_rad_arcsec = round_number(float(eff_rad) * pxsc_zcal_const[0],3)
    eff_rad_arcsec_err = round_number(float(eff_rad_err) * pxsc_zcal_const[0],3)
    
    # Effective radius in kpc
    #moving distances from kpc to arcsec
    
    if z != 0.0:
        zscale = cosmological_scale(z)
    else:
        zlocal = zlocal_correction(galaxy_name)
        zscale = cosmological_scale(zlocal)
        
    
    
    eff_rad_kpc = round_number(eff_rad_arcsec * zscale,3)
    eff_rad_kpc_err = round_number(eff_rad_arcsec_err * zscale,3)

    
    ser_index = hdr['1_N'].split(' ')[0]
    ser_index_err = hdr['1_N'].split(' ')[2]


    if '*' in ser_index:
        ser_index = np.nan
        ser_index_err = np.nan
    
    ax_rat = hdr['1_AR'].split(' ')[0]
    ax_rat_err = hdr['1_AR'].split(' ')[2]

    if '*' in ax_rat:
        ax_rat = np.nan
        ax_rat_err = np.nan
    
    pos_ang = hdr['1_PA'].split(' ')[0]
    pos_ang_err = hdr['1_PA'].split(' ')[2]

    if '*' in pos_ang:
        pos_ang = np.nan
        pos_ang_err = np.nan
    
    chi2 = hdr['CHISQ']
    chi2nu = hdr['CHI2NU']
    
    # adding all the values to a new column of the dataframe
    df.loc[len(df)]=[z,x_center,x_center_err,y_center,y_center_err,mag,mag_err,eff_rad,eff_rad_err,eff_rad_arcsec,eff_rad_arcsec_err,zscale,eff_rad_kpc,eff_rad_kpc_err,ser_index,ser_index_err,ax_rat,ax_rat_err,pos_ang,pos_ang_err,chi2,chi2nu]
    df = df.astype('float64')
    
    return df,y_center,x_center,mag,eff_rad,ser_index,ax_rat,pos_ang

"""

def parse_galfit_param(header_string):
    """
    Parses a Galfit header string ('value +/- error'),
    removes asterisks if present, and converts to float.
    Returns (np.nan, np.nan) if asterisks are found or conversion fails.
    """
    try:
        parts = header_string.split() # Split by spaces: ['value', '+/-', 'error']
        value_str = parts[0]
        error_str = parts[2]

        cleaned_value_str = value_str.strip('*')
        cleaned_error_str = error_str.strip('*')

        value_float = float(cleaned_value_str)
        error_float = float(cleaned_error_str)

        return value_float, error_float
    
    except (IndexError, ValueError, TypeError):
        # If anything fails (unexpected format, cannot convert to float), return NaN
        return np.nan, np.nan

# --- Your modified main function --- #

def galfit_create_dataframe(df, galaxy_name, output_file_path, pxsc_zcal_const):

    print('Creating a dataframe to store results') # Keep simple prints, or consider using logging

    hdr, img, frame_output_name = open_fits(output_file_path, dim=2)

    # Find the redshift in the output file name
    # Ensure z is float for calculations
    try:
        # Use raw string r'...' for regex patterns
        z = float(re.findall(r'z(\d\.\d\d)', output_file_path)[0])
    except (IndexError, ValueError):
        z = np.nan # Handle case if z is not found or not convertible

    # Use the helper function for all parameters
    # Use hdr.get() with a default value in case the key is missing
    x_center, x_center_err = parse_galfit_param(hdr.get('1_XC', '*nan* +/- *nan*'))
    y_center, y_center_err = parse_galfit_param(hdr.get('1_YC', '*nan* +/- *nan*'))
    mag, mag_err = parse_galfit_param(hdr.get('1_MAG', '*nan* +/- *nan*'))
    eff_rad, eff_rad_err = parse_galfit_param(hdr.get('1_RE', '*nan* +/- *nan*'))
    ser_index, ser_index_err = parse_galfit_param(hdr.get('1_N', '*nan* +/- *nan*'))
    ax_rat, ax_rat_err = parse_galfit_param(hdr.get('1_AR', '*nan* +/- *nan*'))
    pos_ang, pos_ang_err = parse_galfit_param(hdr.get('1_PA', '*nan* +/- *nan*'))

    # Calculate derived values (ensure round_number handles NaNs or check before calling)
    if not np.isnan(eff_rad):
        eff_rad_arcsec = round_number(eff_rad * pxsc_zcal_const[0], 3)
        eff_rad_arcsec_err = round_number(eff_rad_err * pxsc_zcal_const[0], 3)

        # Calculate zscale (ensure cosmological_scale/zlocal_correction return float!)
        if z != 0.0 and not np.isnan(z):
             zscale = cosmological_scale(z)
        elif not np.isnan(z): # Implies z == 0.0
             zlocal = zlocal_correction(galaxy_name) # Assumes this returns an effective z > 0
             zscale = cosmological_scale(zlocal)
        else: # If z is NaN
             zscale = np.nan

        if not np.isnan(zscale):
             eff_rad_kpc = round_number(eff_rad_arcsec * zscale, 3)
             eff_rad_kpc_err = round_number(eff_rad_arcsec_err * zscale, 3)
        else:
             eff_rad_kpc = np.nan
             eff_rad_kpc_err = np.nan

    else: # If eff_rad is NaN
        eff_rad_arcsec = np.nan
        eff_rad_arcsec_err = np.nan
        eff_rad_kpc = np.nan
        eff_rad_kpc_err = np.nan
        zscale = np.nan # zscale is only used for eff_rad_kpc

    # Try reading Chi^2, converting to float if possible
    try:
        chi2 = float(hdr.get('CHISQ', np.nan))
        chi2nu = float(hdr.get('CHI2NU', np.nan))
    except (ValueError, TypeError):
        chi2 = np.nan
        chi2nu = np.nan

    # Adding all the values (now mostly floats or NaNs) to the dataframe
    df.loc[len(df)] = [z, x_center, x_center_err, y_center, y_center_err, mag, mag_err,
                       eff_rad, eff_rad_err, eff_rad_arcsec, eff_rad_arcsec_err,
                       zscale, eff_rad_kpc, eff_rad_kpc_err, ser_index, ser_index_err,
                       ax_rat, ax_rat_err, pos_ang, pos_ang_err, chi2, chi2nu]

    return df, y_center, x_center, mag, eff_rad, ser_index, ax_rat, pos_ang


def imfit_init_dataframe():
    
    # creating the Pandas dataframe to store the information
    df = pd.DataFrame(columns = ['z',
                                 'X_cent', 
                                 'X_cent_err', 
                                 'Y_cent',
                                 'Y_cent_err', 
                                 'I_e','I_e_err',
                                 'Eff_rad',
                                 'Eff_rad_err',
                                 'Eff_rad_arcsec',
                                 'Eff_rad_arcsec_err',
                                 'kpc_per_arcsec',
                                 'Eff_rad_kpc',
                                 'Eff_rad_kpc_err',
                                 'Ind','Ind_err',
                                 'Ax_rat',
                                 'Ax_rat_err',
                                 'Pos_ang',
                                 'Pos_ang_err',
                                 'Chi2',
                                 'Chi2nu'])

    return df


def imfit_create_dataframe(df,galaxy_name,redshift,param_file,pxsc_zcal_const):

    datos = {}

    with open(param_file, "r") as f:
        for linea in f:
            # Buscar valores clave como Best-fit, AIC, BIC...
            match = re.search(r"^\#\s+(Best-fit value|Reduced value|AIC|BIC):\s+([\d\.]+)", linea)
            if match:
                datos[match.group(1)] = float(match.group(2))
            
            # Buscar parámetros con valores y errores
            match = re.search(r"^(\S+)\s+([\d\.]+)\s+# \+/- ([\d\.]+)", linea)
            if match:
                datos[match.group(1)] = (float(match.group(2)), float(match.group(3)))  # (valor, error)
                
                
    # find the redshif in the ouput file name
    z = redshift
    
    # output values from galfit and their errors
    x_center = datos['X0'][0]
    x_center_err = datos['X0'][1]
    
    y_center =  datos['Y0'][0]
    y_center_err = datos['Y0'][1]
    
    #mag = hdr['1_MAG'].split(' ')[0]
    #mag_err = hdr['1_MAG'].split(' ')[2]
    
    if 'I_e' in datos.keys():
        I_e = datos['I_e'][0]
        I_e_err = datos['I_e'][1]
    else:
        I_e = np.nan
        I_e_err = np.nan
        
    
    eff_rad = datos['r_e'][0]
    eff_rad_err = datos['r_e'][1]

    #if '*' in eff_rad:
    #    eff_rad = np.nan
    #    eff_rad_err = np.nan
    
    # Effective radius in arcsec
    eff_rad_arcsec = round_number(float(eff_rad) * pxsc_zcal_const[0],3)
    eff_rad_arcsec_err = round_number(float(eff_rad_err) * pxsc_zcal_const[0],3)
    
    # Effective radius in kpc
    #moving distances from kpc to arcsec
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(z)
    
    # Extract the unit to create a Quantity inside the if
    kpc_unit = kpc_per_arcsec.unit

    #eff_rad_kpc = eff_rad_arcsec * kpc_per_arcsec.value    
    #eff_rad_kpc_err = eff_rad_arcsec_err * kpc_per_arcsec.value
    
    #if z == '0.00':
        
    # Obtaining a Quantity
    kpc_per_arcsec = zlocal_correction(galaxy_name)*kpc_unit
    
    eff_rad_kpc = round_number(eff_rad_arcsec * kpc_per_arcsec.value,3)
    eff_rad_kpc_err = round_number(eff_rad_arcsec_err * kpc_per_arcsec.value,3)

    ser_index = datos['n'][0]
    ser_index_err = datos['n'][1]

    #if '*' in ser_index:
    #    ser_index = np.nan
    #    ser_index_err = np.nan
        
    ell = datos['ell'][0]
    ell_err =  datos['ell'][1]
    
    ax_rat = round_number(1-ell,3)
    ax_rat_err = ell_err

    #if '*' in ax_rat:
    #    ax_rat = np.nan
    #    ax_rat_err = np.nan
    
    pos_ang = datos['PA'][0]
    pos_ang_err = datos['PA'][1]

    #if '*' in pos_ang:
    #    pos_ang = np.nan
    #    pos_ang_err = np.nan
    
    chi2 = datos['Best-fit value']
    chi2nu = datos['Reduced value']
    
    # adding all the values to a new column of the dataframe
    df.loc[len(df)]=([z,x_center,x_center_err,y_center,y_center_err,I_e,I_e_err,
                      eff_rad,eff_rad_err,eff_rad_arcsec,eff_rad_arcsec_err,
                      kpc_per_arcsec.value,eff_rad_kpc,eff_rad_kpc_err,
                      ser_index,ser_index_err,
                       ax_rat,ax_rat_err,round_number(pos_ang-180,3),pos_ang_err,
                      chi2,chi2nu])
    df = df.astype('float64')
    
    return df,y_center,x_center,I_e,eff_rad,ser_index,ell,pos_ang

