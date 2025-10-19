import os
import shutil
import subprocess
import pdb

# Moving all galfit scripts created to new folders
# to mantein order within the directories

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

# looping inside all directories and files
for (dirpath, dirnames, filenames) in os.walk(folder_new_version_path):
	
	# loop through all subdirectories
	for subdir in dirnames:
		
		# checking if the directorie is the one we want
		if os.path.isfile(f'{folder_new_version_path}/{subdir}/fit.log') == True:
		
			# creating a folder for galfit. scripts
			galfit_folder = f'{subdir}_galfit_scripts'
			
			# if the folder is created then we remove it to overwrite data
			if os.path.isdir(f'{folder_new_version_path}/{subdir}/{galfit_folder}') == True:
			
				shutil.rmtree(f'{folder_new_version_path}/{subdir}/{galfit_folder}')
			
			# if not is created then create the folder
			os.mkdir(f'{folder_new_version_path}/{subdir}/{galfit_folder}')
			
			# looping throught all files in the directorie
			#pdb.set_trace()
			for file in os.listdir(f'{folder_new_version_path}/{subdir}'):

				# selecting just the galfit scripts
				if 'galfit.' in file or 'fit.' in file or '.script' in file:
					
					# moving all files to the new folder
					os.replace(f'{folder_new_version_path}/{subdir}/{file}', f'{folder_new_version_path}/{subdir}/{galfit_folder}/{file}')
