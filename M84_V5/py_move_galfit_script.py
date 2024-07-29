import os
import shutil


# Moving all galfit scripts created to new folders
# to mantein order within the directories

cwd = os.getcwd()

# looping inside all directories and files
for (dirpath, dirnames, filenames) in os.walk(cwd):
	
	# loop through all subdirectories
	for subdir in dirnames:
		
		# checking if the directorie is the one we want
		if os.path.isfile(f'{cwd}/{subdir}/fit.log') == True:
		
			# creating a folder for galfit. scripts
			galfit_folder = f'{subdir}_galfit_scripts'
			
			# if the folder is created then we remove it to overwrite data
			if os.path.isdir(f'{cwd}/{subdir}/{galfit_folder}') == True:
			
				shutil.rmtree(f'{cwd}/{subdir}/{galfit_folder}')
			
			# if not is created then create the folder
			os.mkdir(f'{cwd}/{subdir}/{galfit_folder}')
			
			# looping throught all files in the directorie
			for file in os.listdir(subdir):
				
				# selecting just the galfit scripts
				if 'galfit.' in file or 'fit.' in file:
					
					# moving all files to the new folder
					os.replace(f'{cwd}/{subdir}/{file}', f'{cwd}/{subdir}/{galfit_folder}/{file}')
