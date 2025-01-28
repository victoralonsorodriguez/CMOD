# This script moves the plots created to a new folder
# inside the main galaxy directory

import os
import pdb
import shutil
import subprocess
from datetime import datetime

# getting the current working directory
cwd = os.getcwd()
previus_cwd = cwd + '/..'
global_cwd = previus_cwd + '/../..'

# Galaxy's name
galaxy = cwd.split('/')[-2]

# Path to the folder with the datacubes
fits_path = f'{previus_cwd}/{galaxy}_original_fits'

# Checking the latest version run
version_file_path = f'{previus_cwd}/{galaxy}_version.txt'
version_file = open(f'{version_file_path}','r')
last_line_version_file = str(subprocess.check_output(['tail', '-1', version_file_path]))[2:-1]

folder_new_version_name = last_line_version_file.split(' # ')[0]
folder_new_version_path = f'{previus_cwd}/{folder_new_version_name}'

# Depending on the folder name a large PSF or a
# small psf will be created
folder_name_elements = folder_new_version_name.split('_')
version = folder_name_elements[1]

# New folder to store the plots
plot_folder_name = f'{galaxy}_{version}_plots'
plot_folder_path = f'{folder_new_version_path}/{plot_folder_name}'

date = datetime.today().strftime('%Y%m%d')

plot_folder_global_name = f'galfit_plots_global'
plot_folder_global_path = f'{global_cwd}/{plot_folder_global_name}'

plot_folder_path_dest = f'{plot_folder_global_path}/{plot_folder_name}'

# Creating the global folder if it is not created yet
if os.path.isdir(plot_folder_global_path) == False:
    
    os.mkdir(plot_folder_global_path)

# Removing the individual plot folder if is already created
if os.path.isdir(plot_folder_path_dest) == True:

    shutil.rmtree(plot_folder_path_dest)

# Moving the plot folder
shutil.copytree(plot_folder_path, plot_folder_path_dest)  
