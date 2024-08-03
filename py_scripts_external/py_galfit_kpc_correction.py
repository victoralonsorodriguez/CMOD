
# This scripts avoids the 0 kpc galaxy radius at z=0
def kpc_correction(galaxy):

    # Dictionry to store the z=0 kpc/'' factor correction
    scale_dict = {'ESO498G05':0.192,
                    'IC719':0.154,
                    'IC2051':0.125,
                    'M84':0.095,
                    'NGC0289':0.096,
                    'NGC307':0.250,
                    'NGC788':0.266,
                    'NGC1309':0.140,
                    'NGC1440':0.105,
                    'NGC1553':0.075,
                    'NGC3393':0.285,
                    'NGC3783':0.226,
                    'NGC4418':0.174,
                    'NGC5806':0.110,
                    'NGC6958':0.176
                    }
                    
    return scale_dict[galaxy]
