
# This script removes the py_scripts_external from
# all the galaxies when they are dividen in 
# Earlt, Late and Other Type

for dir in ../*
do 
    if [ -d $dir ] && [ $dir != ../py_scripts_external ] && [[ $dir != *"plot"* ]]
        then
        cd ./$dir
        for subdir in */
        do
            echo $subdir
            cd ./$subdir
            for subsubdir in */
            do
                if  [ -d $subsubdir ] && [[ $subsubdir != *"_original_fits/" ]]
                then
                    echo Removing $subsubdir
                    rm -r $subsubdir
                    rm *version.txt
                fi
            done
            cd ../
        done
        cd ../py_scripts_external
    fi
done