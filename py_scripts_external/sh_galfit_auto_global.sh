
for argument in "$@"
do
    case $argument in
        "-P")
            psf_type=-P
            ;;
        "-p") 
            psf_type=-p
            ;;
        "-o")
            version=-o
            ;;
        "-m")
            version=-m
            ;;
        "-mi")
            version=-mi
            ;;      
        "-i")
            version=-i
            ;;  
        *)
            psf_type=-P
            version=-mi
    esac
done 

if [ -z ${psf_type+x} ]
then 
    psf_type=-P
fi

if [ -z ${version+x} ]
then 
    version=-mi
fi

for dir in $ *
do
    if [ -d $dir ] && [ $dir != py_scripts_external ]
    then 
        echo $'\n'Starting Galfit analysis for galaxy $dir
        cd ./$dir/py_scripts_external
        ./sh_galfit_auto.sh $psf_type $version
        echo $'\n'Full Galfit analysis is finished for galaxy $dir $'\n'
        cd ../.. 
    fi
done