

from py_round_number import round_number


# Changing from pixeles to kpc 
def px_to_kpc(px,inst):

    arcsec = px * inst 
    kpc = arcsec_to_kpc(arcsec)

    return round_number(kpc,3)

def kpc_to_px(kpc,inst):
    
    arcsec = kpc_to_arcsec(kpc)
    px = arcsec / inst 
    
    return round_number(px,3)

def arcsec_to_kpc(arcsec):
    
    kpc = arcsec * 0.12
    
    return round_number(kpc,3)

def kpc_to_arcsec(kpc):
    
    arcsec = kpc / (0.12)
    
    return round_number(arcsec,3)