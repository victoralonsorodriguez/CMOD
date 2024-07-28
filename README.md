# CMOD - Chromatic Surface Brightness MODulation

Ths repository contains the codes created to study the CMOD effect from a paractical point of view.

The main branch contains the final files of every development branch correponding to each new version. 

## Version change log

### Version 2.03 - unite_version branch
This version simplify the files required to execute the different versions of the code. With this version all the diferent modes of analysis are contained in the same file 'py_galfit_3D.py'. To execute the modes it is neccessary to indicate them as anrgument in the comment line. For example, to run the original mode it should be as 'python3 py_galfit_3D.py original'. If no argument is given, the default mode is the 'medfilt_version' since is the most accurated mode. 

### Version 2.02 - plots_from_csv_folder
Now the plots are created from the csv folder for each galaxy. New plots of the magnitude along with the Sérsic Index and Effective Radiud are created. Also the center position along with the Axis Ratio and the Position Angle.

### Version 2.01 - csv_files branch
This version includes a new folder to store the csv files created during the galfit analysis. This folder is included inside the main galaxy directory names as 'galaxyname_csv'.

## Version 1.00
This version was presented as the final code for the BSc Thesis. It allows to analyse the galaxies with galfit and create some plots.

