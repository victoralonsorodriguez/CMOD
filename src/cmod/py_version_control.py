import pdb

import os
import re
import shutil
import numpy as np
from datetime import datetime

from py_create_folder import create_folder

def version_directory(cwd,base_name,original_dir='.',create_dir=False):

    pattern = re.compile(rf'{re.escape(base_name)}_V(\d+)')
    versions = []

    # Searching for version
    for dir_name in os.listdir(cwd):
        match = re.match(pattern, dir_name)
        if match:
            versions.append(int(match.group(1)))
    
    # If no version is found
    if len(versions)==0:
        versions.append(0)
    
    # Determining the new version
    versions = np.array(versions)
    previus_version = max(np.nanmax(versions),0)
    next_version = previus_version + 1
    
    # Creating the folder
    next_dir_path = f'{original_dir}_V{next_version}'
    if create_dir == True:
        create_folder(next_dir_path)
    
    return int(previus_version)


def version_file(galaxy_name, version, analysis_version, psf_size,archivo,analysis_programe):

    datetime_current = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    date, time = datetime_current.split(' ')
    
    nueva_linea = (f'VERSION={version} PROGRAMME={analysis_programe} DATE={date} TIME={time} ' 
                   f'ANALYSIS_MODE={analysis_version} PSF={psf_size}\n')

    # Comprobar si el archivo ya existe
    existe = os.path.exists(archivo)

    with open(archivo, 'a' if existe else 'w') as f:
        # Si el archivo no existe, escribir el encabezado
        if not existe:
            f.write(f'#-------GALAXY {galaxy_name}-------#\n\n')
        f.write(nueva_linea)


def version_file_last(file):

    if not os.path.exists(file):   
        return 0

    versions = []
    with open(file, 'r') as f:
        for line in f:
            match = re.search(r'VERSION=(\d+)', line)
            if match:
                versions.append(int(match.group(1)))
    
    max_version = max(versions)

    return int(max_version)