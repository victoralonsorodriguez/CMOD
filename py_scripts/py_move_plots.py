
# This script moves the plots created to a new folder
# inside the main galaxy directory

import os
import shutil


cwd = os.getcwd()
galaxy = (cwd.split('/')[-1])

# New folder to store the plots
plot_folder = f'{galaxy}_plots'

if os.path.isdir(plot_folder) == True:
    
    shutil.rmtree(f'{cwd}/{plot_folder}')

os.mkdir(plot_folder)

# Moving the plots
for file in os.listdir(cwd):

    if '.pdf' in file and 'plot' in file and ('galfit' in file or 'ratios' in file):
            
            file_original_path = os.path.join(cwd, file)
            file_new_path = os.path.join(cwd,f'{plot_folder}',file)
            
            shutil.copyfile(file_original_path,file_new_path)
            os.remove(file_original_path)            