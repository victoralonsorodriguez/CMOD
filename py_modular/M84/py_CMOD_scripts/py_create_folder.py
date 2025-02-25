
import os
import shutil

def create_folder(folder_path,
                  overwrite = False):
    
    # Check if the folder already exists
    if os.path.isdir(folder_path) == True:
        
        # If overwrite is selected
        if overwrite == True:
            shutil.rmtree(folder_path)
            os.mkdir(folder_path)

        # If not, pass
        else:
            pass
    
    # If the folder is not already created
    else:
        os.mkdir(folder_path)
        