
# This script will create a txt file
# in which it is stored the vervsion of the analysis

import os
import shutil
import re
import subprocess

import pdb



cwd = os.getcwd()
previus_cwd = cwd + '/..'

galaxy = cwd.split('/')[-2]

version_folder_list = []
last_run_version = 0

version_file_path = f'{previus_cwd}/{galaxy}_version_file.txt'

if os.path.isfile(version_file_path) == False:
    
    version_file = open(f'{version_file_path}','a+')
    version_file.write(f'This file contains the run versions for the galaxy {galaxy}\n\n')
    version_file.write(f'#----------VERSIONS RUN----------#\n\n')

version_file = open(f'{version_file_path}','a+')



last_line_version_file = str(subprocess.check_output(['tail', '-1', version_file_path]))[2:-1]

last_version_file_name = last_line_version_file.split(' # ')[0]

if len(re.findall(f'{galaxy}_V'+r'\d',str(last_version_file_name))) != 0:

    last_version_file_num = last_version_file_name.split('_')[-1].split('V')[1]

    last_run_version = last_version_file_num

    print(f'Last registered version was V{last_run_version}')

else:

    last_version_file_name = f'{galaxy}_V0'
    last_version_file_num = last_version_file_name.split('_')[-1].split('V')[1]
    last_run_version = last_version_file_num

    print('No registered versions on the file')



for folder in sorted(os.listdir(previus_cwd)):

    if len(re.findall(f'{galaxy}_V'+r'\d',str(folder))) != 0:

        version_folder_list.append(list)
        folder_version_num = folder.split('_')[-1].split('V')[1]
        
        if int(folder_version_num)>int(last_run_version):

            last_run_version = folder_version_num
            last_version_file_name = folder


folder_new_version_name = last_version_file_name.split(f'{last_run_version}')[0]+f'{int(last_run_version)+1}'

os.mkdir(f'{previus_cwd}/{folder_new_version_name}')

version_file.write(f'{folder_new_version_name} # \n')