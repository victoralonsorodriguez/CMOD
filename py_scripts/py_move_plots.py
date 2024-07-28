import os
import shutil


# This script moves the plots created to a new folder
# inside the main galaxy directory


cwd = os.getcwd()
galaxy = (cwd.split('/')[-1])

# For plots with the data
plot_folder = f'{galaxy}_plots_galfit'

if os.path.isdir(plot_folder) == True:
    
    shutil.rmtree(f'{cwd}/{plot_folder}')

os.mkdir(plot_folder)

for file in os.listdir(cwd):

    if '.pdf' in file and 'multiplot' in file and 'galfit' in file:
            
            file_original_path = os.path.join(cwd, file)
            file_new_path = os.path.join(cwd,f'{plot_folder}',file)
            
            shutil.copyfile(file_original_path,file_new_path)
            os.remove(file_original_path)            

# For plots with ratios
plot_folder = f'{galaxy}_plots_ratios'

if os.path.isdir(plot_folder) == True:
    
    shutil.rmtree(f'{cwd}/{plot_folder}')

os.mkdir(plot_folder)

for file in os.listdir(cwd):

    if '.pdf' in file and 'multiplot' in file and 'ratios' in file:
            
            file_original_path = os.path.join(cwd, file)
            file_new_path = os.path.join(cwd,f'{plot_folder}',file)
            
            shutil.copyfile(file_original_path,file_new_path)
            os.remove(file_original_path)
