

import numpy as np

from py_round_number import round_number


# Changing from degrees to radian
def deg_to_rad(deg):
    
    rad = (deg) * (np.pi / 180)

    return round_number(rad,3)

def rad_to_deg(rad):
    
    deg = (rad) * (180 / np.pi)
        
    return round_number(deg,3)


def rad_to_deg_abs(rad):
    
    deg = (rad) * (180 / np.pi)
        
    return round_number(deg,3)