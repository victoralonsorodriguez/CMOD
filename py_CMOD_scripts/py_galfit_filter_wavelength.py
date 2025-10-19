def filter_wl_dict():
    
    # The names of the filter systems in alphabetic order
    filter_system_names = ['Euclid','JWST','Hubble','Johnson_UVRIJHK']
    filter_system_corrected_names = ['Euclid','JWST','HST','SDSS']
    
    # Dictionary with the corresponding wavelength of each filter
    # http://svo2.cab.inta-csic.es/theory/fps/
    filter_wavelenght_dict = {'EucHab':14971.70, 'EucJab':11522.58 , 'EucVISab':4970.77,
                'EucYab':9381.52,
                'F070W':7039.12, 'F090W':9021.53, 'F115W':11542.61, 
                'F140M':14053.23, 'F150W':15007.44, 'F162M':16272.47, 
                'F182M':18451.67, 'F200W':19886.48, 'F210M':20954.51, 
                'F277W':27617.40, 'F356W':35683.62, 'F444W':44043.15,
                'Hve':16396.38, 'Ive':8657.44, 'Jve':12317.30, 'Kve':21735.85,
                'Rve':6819.05, 'Uve':3511.89, 'Vve':5501.40,
                'HSTF300W':2984.47, 'HSTF435W':4323.09, 'HSTF450W':4556.22,
                'HSTF475W':4775.73, 'HSTF555W':5355.74, 'HSTF569W':5644.08, 
                'HSTF606W':5887.08, 'HSTF791W':7875.56, 
                'HSTF814W':8039.03}
    
    # Dictionary with the correct name to plot for each filter
    filter_names_dict = {'EucHab':'HE', 'EucJab':'JE' , 'EucVISab':'IE',
                'EucYab':'YE',
                'F070W':'F070W', 'F090W':'F090W', 'F115W':'F115W', 
                'F140M':'F140M', 'F150W':'F150W', 'F162M':'F162M', 
                'F182M':'F182M', 'F200W':'F200W', 'F210M':'F210M', 
                'F277W':'F277W', 'F356W':'F356W', 'F444W':'F444W',
                'Hve':'H', 'Ive':'I', 'Jve':'J', 'Kve':'K',
                'Rve':'R', 'Uve':'U', 'Vve':'V',
                'HSTF300W':'F300W', 'HSTF435W':'F435W', 'HSTF450W':'F450W',
                'HSTF475W':'F475W', 'HSTF555W':'F555W', 'HSTF569W':'F569W', 
                'HSTF606W':'F606W', 'HSTF791W':'F791W', 
                'HSTF814W':'F814W'}
                
                
    return filter_system_corrected_names,filter_wavelenght_dict,filter_names_dict
