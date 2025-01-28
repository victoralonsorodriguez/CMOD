
# This script moves the plots created to a new folder
# inside the main galaxy directory

import os
import pdb
import shutil
import subprocess

# getting the current working directory
cwd = os.getcwd()
previus_cwd = cwd + '/..'

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

if os.path.isdir(plot_folder_path) == True:
    
    shutil.rmtree(f'{plot_folder_path}')

os.mkdir(plot_folder_path)

# Moving the plots
for file in os.listdir(folder_new_version_path):

    if '.pdf' in file and 'plot' in file and ('galfit' in file or 'ratios' in file):
            
            file_original_path = os.path.join(folder_new_version_path, file)
            file_new_path = os.path.join(f'{plot_folder_path}',file)
            
            shutil.copyfile(file_original_path,file_new_path)
            os.remove(file_original_path)