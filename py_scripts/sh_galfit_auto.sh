
# This shell script executes the entire pipeline for
# the galfit analisys of each galaxy

# Creating the large psf
python3 py_psf.py

# Version to run: original, medfilt
# inout or medfilt_ionut (default)
version=$1
if [ -n "$version" ]
then 
    python3 py_galfit_3D.py $version
else
    python3 py_galfit_3D.py medfilt_inout
fi

# Moving generated files to each filter directoy
python3 py_move_galfit_script.py

# Creating the plots
python3 py_plot_galfit_all.py

# Moving the plots to an specific folder
python3 py_move_plots.py

echo Galfit full analysis is finished with medfilt_inout methods
