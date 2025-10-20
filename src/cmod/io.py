

import sys
import argparse
import re
import ast
import pdb
import os


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
    hdu = fits.open(fits_path)
    hdr = hdu[dim].header
    img = hdu[dim].data
    
    # Returning the header, the image and the name
    return (hdr,img,fits_name)


def load_configuration(conf_file):
    conf_var = {}
    with open(conf_file, "r") as file:
        for line in file:
            line = line.split("#")[0].strip()  # Eliminar comentarios y espacios
            if line:  # Ignorar líneas vacías
                match = re.match(r"^([A-Z_]+)\s*=\s*(.+)", line)  # Buscar VARIABLE = VALOR
                if match:
                    key, value = match.groups()
                    value = value.strip()
                    
                    # Convertir valores: int si es número, lista o tupla si tiene formato correcto
                    try:
                        value = ast.literal_eval(value)  # Convierte '0' a 0, '["x"]' a lista, etc.
                    except (ValueError, SyntaxError):
                        pass  # Si no es un número o estructura reconocible, se deja como string
                    
                    conf_var[key] = value
    
    return conf_var

def argparse_values(phase):
    
    # All possible flags
    parser = argparse.ArgumentParser()
    parser.add_argument('-conf','--confifile', type=str, 
                        help=(f'Path to the configuration file. If this flag is marked,\n'
                              f'configuration file will be use as deafult option'))
    
    # Analysis configuration options
    parser.add_argument('-ap','--analysis_programme', type=int, nargs='?', const=0, default=0,
                         help=(f'Programme to compute the analysis:\n'
                              f'=0: Galfit (Default);\n'
                              f'=1: Imfit'))
    parser.add_argument('-am','--analysis_mode', type=int, nargs='?', const=3, default=3,
                         help=(f'Version to compute the analysis:\n'
                              f'=0: Original;\n'
                              f'=1: Medfilt;\n'
                              f'=2: Inout;\n'
                              f'=3: Medinout (Default)'))
    parser.add_argument('-ar','--analysis_range', type=int, nargs='?', const=0, default=0,
                         help=(f'Until which frame is analyzed:\n'
                              f'=0: Full range (Default);\n'
                              f'=N: Until frame N'))
    parser.add_argument('-psf','--psf', type=int, nargs='?', const=0, default=0,
                         help=(f'Point Spread Function:\n'
                              f'=0: Auto Large (Default);\n'
                              f'=1: Auto Small;\n'
                              f'=/path/: Specific'))
    parser.add_argument('-ic','--initial_conditions', type=int, nargs='?', const=0, default=0,
                         help=(f'Initial conditions determination:\n'
                              f'=0: Auto (Default);\n'
                              f'=/path/: Specific IC file path'))
    
    # Version control option
    parser.add_argument('-vconti','--version_continuation', type=int, nargs='?', const=1, default=1,
                         help=(f'Continue analisys if interrupted:\n'
                              f'=0: No;\n'
                              f'=1: Yes (Default)'))   
    parser.add_argument('-vcontr','--version_control', type=int, nargs='?', const=1, default=1,
                         help=(f'Keep track of the runned version:\n'
                              f'=0: No;\n'
                              f'=1: Yes (Default)'))   
    
    # Values required
    parser.add_argument('-pxsc','--pixel_scale', type=int, nargs='?', const=0.2, default=0.2,
                         help=(f'Intrument pixel scale in arsec/pixel:\n'
                              f'=0.2 (Default)\n'))   
    parser.add_argument('-zcal','--zero_calibration', type=int, nargs='?', const=25, default=25,
                         help=(f'Instrument calibration constant:\n'
                              f'=25 (Default)'))
    parser.add_argument('-cropf','--crop_factor', type=int, nargs='?', const=10, default=10,
                         help=(f'Galaxy center cropping factor:\n'
                              f'=10 (Default)'))         
    
    # Plot configuration options
    parser.add_argument('-ps','--plot_system', type=int, nargs='?', const=0, default=0,
                        help=(f'Systems of filters to plot:\n'
                              f'=0: Automatic, all vailable systems are plotted (Default);\n'
                              f'=["","",...]: list of the specific systems to plot'))
    parser.add_argument('-pf','--plot_filter', type=int, nargs='?', const=0, default=0,
                        help=(f'Filters to plot:\n'
                              f'=0: Automatic, all vailable filters are plotted (Default);\n'
                              f'=["","",...]: list of the specific filters to plot'))
    parser.add_argument('-pmag','--plot_magnitude', type=int, nargs='?', const=0, default=0,
                        help=(f'Magnitudes to plot:\n'
                              f'=0: Automatic, all vailable magnitudes are plotted (Default);\n'
                              f'=[("",""),...]: list of tuples of the specific magnitudes to plot'))
    parser.add_argument('-pm','--plot_mode', type=int, nargs='?', const=2,default=2,
                        help=(f'Mode of plot data:\n'
                              f'=0: Filter, each filter will be plotted for a pair of magnitudes in the same sheet;\n'
                              f'=1: Magnitude, all the magnitudes will be plot for each filter system;\n'
                              f'=2: Both modes are plotted (Default);\n'))
    parser.add_argument('-psm','--plot_save_mode', type=int, nargs='?', const=0, default=0,
                        help=(f'Mode of saving plot data:\n'
                              f'=0: All filters plots as a group (Default);\n'
                              f'=1: All filters individually\n'))
    parser.add_argument('-pmc','--plot_max_columns', type=int, nargs='?', const=12, default=12,
                        help=(f'Maximun number of columns for each filter:\n'
                              f'=N: Specific;\n'
                              f'=8: (Default)\n'))
    parser.add_argument('-pav','--plot_analysis_version', type=int, nargs='?', const=99, default=99,
                        help=(f'Maximun number of columns for each filter:\n'
                              f'=N: Specific;\n'
                              f'=latest: (Default)\n'))
    
    # Obtaining the used flags
    pargs = parser.parse_args()
    sargs = sys.argv
    
    # If there is a configuration file then it is used
    # as default configurator
    if pargs.confifile != None:
        
        # Load configuration file
        conf_var = load_configuration(pargs.confifile)
        
        # Analysis configuration options
        analysis_programme = conf_var['ANALYSIS_PROGRAMME']
        analysis_mode = conf_var['ANALYSIS_MODE']
        psf = conf_var['PSF']
        initial_conditions = conf_var['INITIAL_CONDITIONS']
        
        if '-ap' in sargs or '--analysis_programme' in sargs:
            analysis_programme = pargs.analysis_programme
        if '-am' in sargs or '--analysis_mode' in sargs:
            analysis_mode = pargs.analysis_mode
        if '-ar' in sargs or '--analysis_range' in sargs:
            analysis_range = pargs.analysis_range
        if '-psf' in sargs or '--psf' in sargs:
            psf = pargs.psf
        if '-ic' in sargs or '--initial_conditions' in sargs:
            initial_conditions = pargs.initial_conditions    
    
        # Version control option
        version_continuation = conf_var['VERSION_CONTINUATION']
        version_control = conf_var['VERSION_CONTROL']
        
        if '-vconti' in sargs or '--version_continuation' in sargs:
            version_continuation = pargs.version_continuation
        if '-vcontr' in sargs or '--version_control' in sargs:
            version_control = pargs.version_control       
            
        # Values required
        pixel_scale = conf_var['PIXEL_SCALE']
        zero_calibration = conf_var['ZERO_CALIBRATION_CONSTANT']
        crop_factor = conf_var['CROP_FACTOR']
        
        if '-pxsc' in sargs or '--pixel_scale' in sargs:
            pixel_scale = pargs.pixel_scale
        if '-zcal' in sargs or '--zero_calibration' in sargs:
            zero_calibration = pargs.zero_calibration   
        if '-cropf' in sargs or '--crop_factor' in sargs:
            crop_factor = pargs.crop_factor           
        
        
        # Plot configuration options
        system_list_mode = conf_var['PLOT_SYSTEM']
        filters_list_mode = conf_var['PLOT_FILTER']
        mag_pairs_mode = conf_var['PLOT_MAG']
        plot_mode = conf_var['PLOT_MODE']
        plot_save_mode = conf_var['PLOT_SAVE_MODE']
        plot_max_columns = conf_var['PLOT_MAX_COLUMNS']
        plot_analysis_version = conf_var['PLOT_ANALYSIS_VERSION']
        
        # If one argumen is parsed, then overwritte 
        # config file value
        if '-ps' in sargs or '--plot_system' in sargs:
            system_list_mode = pargs.plot_system
        if '-pf' in sargs or '--plot_filter' in sargs:
            filters_list_mode = pargs.plot_filter
        if '-pmag' in sargs or '--plot_magnitude' in sargs:
            mag_pairs_mode = pargs.plot_magnitude
        if '-pm' in sargs or '--plot_mode' in sargs:
            plot_mode = pargs.plot_mode
        if '-psm' in sargs or '--plot_save_mode' in sargs:
            plot_save_mode = pargs.plot_save_mode
        if '-pmc' in sargs or '--plot_max_columns' in sargs:
            plot_max_columns = pargs.plot_max_columns
        if '-pav' in sargs or '--plot_analysis_version' in sargs:
            plot_analysis_version = pargs.plot_analysis_version
    
    # If no config file is used, then default values
    # values or the indicated one are used
    else:
        
        # Analysis configuration options
        analysis_programme = pargs.analysis_programme
        analysis_mode = pargs.analysis_mode
        analysis_range = pargs.analysis_range
        psf = pargs.psf
        initial_conditions = pargs.initial_conditions    
        
        # Version control option
        version_continuation = pargs.version_continuation
        version_control = pargs.version_control    
        
        # Values required
        pixel_scale = pargs.pixel_scale
        zero_calibration = pargs.zero_calibration   
        crop_factor = pargs.crop_factor   
        
        # Plot configuration options
        system_list_mode = pargs.plot_system
        filters_list_mode = pargs.plot_filter
        mag_pairs_mode = pargs.plot_magnitude
        plot_mode = pargs.plot_mode
        plot_save_mode = pargs.plot_save_mode
        plot_max_columns = pargs.plot_max_columns
        plot_analysis_version = pargs.plot_analysis_version
        
    
    if analysis_programme == 0:
        analysis_programme = 'Galfit'
    else:
        analysis_programme = 'Imfit'
        
    if analysis_mode == 0:
        analysis_mode = 'original'
    elif analysis_mode == 1:
        analysis_mode = 'medfilt'
    elif analysis_mode == 2:
        analysis_mode = 'inout'
    elif analysis_mode == 3:
        analysis_mode = 'medinout'
        
    if psf == 0:
        psf = 'large'
    elif psf == 1:
        psf = 'small'
    elif psf == 'str':
        psf = 'path'
        
    if plot_analysis_version == 99:
        plot_analysis_version = 'latest'

    if phase == 'analysis':
        return (analysis_programme,analysis_mode,analysis_range,
                version_continuation,version_control,psf,initial_conditions,
                pixel_scale,zero_calibration,crop_factor)
        
    elif phase == 'plot':
        return (analysis_programme,system_list_mode,filters_list_mode,
                mag_pairs_mode,plot_mode,plot_save_mode,plot_max_columns,
                plot_analysis_version)
        
        
        
def extract_filter_name(fits_filename):
    """
    Extracts the filter name from a CMOD FITS filename.

    Assumes the filter name is the last part of the filename
    before the '.fits' extension, separated by underscores.
    Example: 'n1309_VBIN020_SL_zSimJ_EucHab.fits' -> 'EucHab'

    Args:
        fits_filename (str): The full filename (e.g., 'galaxy_..._filter.fits').

    Returns:
        str: The extracted filter name, or None if the format is unexpected.
    """
    try:
        # Get the base name without the extension
        base_name = os.path.splitext(fits_filename)[0]
        # Split by underscore and take the last part
        filter_name = base_name.split('_')[-1]
        return filter_name
    except IndexError:
        # Handle cases where splitting might fail (e.g., no underscores)
        print(f"Warning: Could not extract filter name from '{fits_filename}'. Unexpected format.")
        return None