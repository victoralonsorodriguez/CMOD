
# This script removes the py_scripts_external from
# all the galaxies when they are dividen in 
# Earlt, Late and Other Type

for dir in ../*
do 
    if [ -d $dir ] && [ $dir != ../py_scripts_external ] && [[ $dir != *"plot"* ]]
        then
        echo $dir
        cd ./$dir
        rm sh_galfit_auto_global.sh
        for subdir in */
        do
            echo $subdir
            rm -r ./$subdir/py_scripts_external
        done
        cd ../py_scripts_external
    fi
done