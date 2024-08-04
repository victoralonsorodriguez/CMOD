
# This script will create a txt file
# in which it is stored the vervsion of the analysis

import os
import sys
import shutil
import re
import subprocess
from datetime import datetime

import pdb



cwd = os.getcwd()
previus_cwd = cwd + '/..'

galaxy = cwd.split('/')[-2]

version_folder_list = []
last_run_version = 0

version_file_path = f'{previus_cwd}/{galaxy}_version.txt'

if os.path.isfile(version_file_path) == False:
    
    version_file = open(f'{version_file_path}','a+')
    version_file.write(f'This file contains the run versions for the galaxy {galaxy}\n\n')
    version_file.write(f'#----------VERSIONS RUN----------#\n\n')

version_file = open(f'{version_file_path}','a+')

# Default version
if len(sys.argv) == 1:

    # By default the running version is the 
    # medfilt_inout
    sys.argv.append('PSF')
    sys.argv.append('medfilt_inout')


# Checking the lastest version from the file
last_line_version_file = str(subprocess.check_output(['tail', '-1', version_file_path]))[2:-1]

last_version_file_name = last_line_version_file.split(' # ')[0]

# Check if the previus version uses a PSF
if 'psf' in str.lower(last_version_file_name).split('_'):

    psf = last_version_file_name.split('_')[-1]

    if len(re.findall(f'{galaxy}_V'+r'\d+'+f'_{psf}',str(last_version_file_name))) != 0:

        last_version_file_num = last_version_file_name.split('PSF')[0].split('V')[-1].split('_')[0]

        last_run_version = last_version_file_num

        #print(f'Last registered version was V{last_run_version}')

    else:

        last_version_file_name = f'{galaxy}_V0'
        last_version_file_num = last_version_file_name.split('_')[-1].split('V')[1]
        last_run_version = last_version_file_num

else:

    if len(re.findall(f'{galaxy}_V'+r'\d+',str(last_version_file_name))) != 0:

        last_version_file_num = last_version_file_name.split('_')[-1].split('V')[1]
        last_run_version = last_version_file_num

    else:

        last_version_file_name = f'{galaxy}_V0'
        last_version_file_num = last_version_file_name.split('_')[-1].split('V')[1]
        last_run_version = last_version_file_num


# Check the last version from the folders
for folder in sorted(os.listdir(previus_cwd)):

    # Check if the folder constains a PSF
    if 'psf' in str.lower(folder).split('_'):

        psf = folder.split('_')[-1]
    
        if len(re.findall(f'{galaxy}_V'+r'\d+'+f'_{psf}',str(folder))) != 0:

            version_folder_list.append(folder)

            folder_version_num = folder.split('PSF')[0].split('V')[-1].split('_')[0]

            last_run_version = last_version_file_num

            if int(folder_version_num)>int(last_run_version):

                last_run_version = folder_version_num
                last_version_file_name = folder

    else: 

        if len(re.findall(f'{galaxy}_V'+r'\d+',str(folder))) != 0:

            version_folder_list.append(folder)

            folder_version_num = folder.split('_')[-1].split('V')[1]
            
            if int(folder_version_num)>int(last_run_version):

                last_run_version = folder_version_num
                last_version_file_name = folder


# Check if the folder exists and is empty then this will be the new version to run
if os.path.isdir(f'{previus_cwd}/{last_version_file_name}') == True and len(os.listdir(f'{previus_cwd}/{last_version_file_name}')) == 0:

    folder_new_version_name = last_version_file_name

# if it is not empty or do not exists then we need to create the new version folder
else:

    # if PSF is given by argument or is an empty argument a PSF will be used
    if sys.argv[1] == 'PSF':

        folder_new_version_name = f'{galaxy}_V{int(last_run_version)+1}'
        os.mkdir(f'{previus_cwd}/{folder_new_version_name}')

    # if psf is given as argument then a psf folder will be created
    else:

        folder_new_version_name = f'{galaxy}_V{int(last_run_version)+1}_{sys.argv[1]}'
        os.mkdir(f'{previus_cwd}/{folder_new_version_name}')


version_file.write(f'\n{folder_new_version_name} # {sys.argv[2]} / {datetime.now()}\n')
