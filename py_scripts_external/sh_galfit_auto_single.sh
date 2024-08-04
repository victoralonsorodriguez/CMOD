
# This shell script executes the entire pipeline for
# the galfit analisys of each galaxy

# Depending on de arguments one action or another
# will be executed
for argument in "$@"
do
    case $argument in
        "-P")
            psf_type=PSF
            ;;
        "-p") 
            psf_type=psf
            ;;
        "-o")
            version=original
            ;;
        "-m")
            version=medfilt
            ;;
        "-mi")
            version=medfilt_inout
            ;;      
        "-i")
            version=inout
            ;;  
        *)
            psf_type=PSF
            version=medfilt_inout
    esac
done 

if [ -z ${psf_type+x} ]
then 
    psf_type=PSF
fi

if [ -z ${version+x} ]
then 
    version=medfilt_inout
fi

echo PSF type: $psf_type // Galfit Version: $version $'\n'
python3 py_version_control.py $psf_type $version

# Creating the large psf
python3 py_psf.py

# Version to run: original, medfilt
# inout or medfilt_ionut (default)
 
python3 py_galfit_3D.py $version

# Moving generated files to each filter directoy
python3 py_move_galfit_script.py

# Creating the plots
python3 py_plot_galfit_all.py

# Moving the plots to an specific folder
python3 py_move_plots.py

echo Galfit full analysis is finished with medfilt_inout methods
