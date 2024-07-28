
# Shell scritp code to automatize the
# full analysis proccess for the
# median filter method

python3 py_psf.py
python3 py_galfit_3D_medfilt.py 
python3 py_move_galfit_script.py
python3 py_plot_galfit_pos_ser_effrad.py
python3 py_plot_ratios_ser_effrad.py
python3 py_move_plots.py
echo Galfit full analysis is finished with medfilt method
