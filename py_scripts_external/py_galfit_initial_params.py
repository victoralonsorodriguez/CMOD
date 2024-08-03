
# Script to store Galfit initial parameters
def initial_params(galaxy):

    if galaxy == 'ESO4958G05':
    
        eff_rad = 73
        ser_ind = 4.0
        ax_rat = 0.85
        pos_ang = -23   

    elif galaxy == 'IC719':
    
        eff_rad = 43
        ser_ind = 1.7
        ax_rat = 0.26
        pos_ang = 52

    elif galaxy == 'IC2051':
    
        eff_rad = 63
        ser_ind = 5.0
        ax_rat = 0.65
        pos_ang = 65

    elif galaxy == 'M84':
    
        eff_rad = 60
        ser_ind = 3.0
        ax_rat = 0.81
        pos_ang = -52

    elif galaxy == 'NGC0289':
    
        eff_rad = 76
        ser_ind = 3.0
        ax_rat = 0.60
        pos_ang = -60

    elif galaxy == 'NGC307':
    
        eff_rad = 30
        ser_ind = 2.7
        ax_rat = 0.51
        pos_ang = 80

    elif galaxy == 'NGC788':
    
        eff_rad = 55
        ser_ind = 2.0
        ax_rat = 0.52
        pos_ang = -75
 
    elif galaxy == 'NGC1309':
    
        eff_rad = 67
        ser_ind = 3.5
        ax_rat = 0.89
        pos_ang = 3

    elif galaxy == 'NGC1440':
    
        eff_rad = 47
        ser_ind = 4.0
        ax_rat = 0.8
        pos_ang = 45

    elif galaxy == 'NGC1553':
    
        eff_rad = 75
        ser_ind = 4.5
        ax_rat = 0.66
        pos_ang = -26

    elif galaxy == 'NGC3393':
    
        eff_rad = 75
        ser_ind = 4.2
        ax_rat = 0.66
        pos_ang = -26

    elif galaxy == 'NGC3783':
    
        eff_rad = 67
        ser_ind = 4.0
        ax_rat = 0.8
        pos_ang = -25

    elif galaxy == 'NGC4418':
    
        eff_rad = 50
        ser_ind = 3.8
        ax_rat = 0.60
        pos_ang = 58

    elif galaxy == 'NGC5806':
    
        eff_rad = 66
        ser_ind = 5.0
        ax_rat = 0.56
        pos_ang = -10

    elif galaxy == 'NGC6958':
    
        eff_rad = 66
        ser_ind = 3.9
        ax_rat = 0.56
        pos_ang = -10

    

    return eff_rad,ser_ind,ax_rat,pos_ang
