
import numpy as np

def round_number(num,dec):
    
    if isinstance(num, (int, float, np.float32, np.float64)):
        num = round(num,dec)
        
    elif isinstance(num, np.ndarray):
        num = np.round(num,dec)   
         
    return num