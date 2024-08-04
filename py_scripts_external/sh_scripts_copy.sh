
# This script copies the py_scripts_external to
# all the galaxies when they are dividen in 
# Earlt, Late and Other Type

for dir in ../*
do 
    if [ -d $dir ] && [ $dir != ../py_scripts_external ]
        then
        echo $dir
        cd ./$dir
        cp ../py_scripts_external/sh_galfit_auto_global.sh .
        for subdir in *
        do
            if [ -d $subdir ] && [ $subdir != py_scripts_external ]
            then
                echo $subdir
                cp -r ../py_scripts_external ./$subdir
            fi
        done
        cd ../py_scripts_external
    fi
done