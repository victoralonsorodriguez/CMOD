

import subprocess
import sys
import pdb

from cmod.utils import Cronometro


if __name__ == '__main__':
    
    arg_flags_list = sys.argv[1:]
    
    cronometro = Cronometro()
    cronometro.iniciar()
    
    #######################################################
    ###------------------CMOD ANALYSIS------------------###
    #######################################################
    
    analysis_flags_list = ['python3','py_analysis.py']+arg_flags_list

    try:
        subprocess.run(analysis_flags_list)
        
    except Exception as e:
        print(f'Error running analysis programe:\n{e}') 
    
    
    ##########################################################
    ###------------------PLOTTING RESULTS------------------###
    ##########################################################
    
    plots_flags_list = ['python3','py_plot_analysis.py']+arg_flags_list
    
    try:
        subprocess.run(plots_flags_list)
        
    except Exception as e:
        print(f'Error running plotting programe:\n{e}')
        
        