# En src/cmod/resampling.py
import numpy as np
from scipy.ndimage import zoom
from astropy.io import fits # Necesario para guardar (por ahora)

def geometric_resample(input_data, factor, order=1):
    """
    Realiza un resampleo geométrico simple de una imagen.

    Args:
        input_data (np.ndarray): Imagen de entrada.
        factor (float): Factor de resampleo. >1 para agrandar, <1 para reducir.
        order (int): Orden de interpolación (0=vecino más cercano, 1=bilineal, 3=cúbica).

    Returns:
        np.ndarray: Imagen resampleada.
    """
    # Asegurarse de que factor_resize sea < 1 para reducir si factor > 1
    # O > 1 para agrandar si factor < 1. 
    # Si factor es "aumentos", el factor de zoom es 1/aumentos.
    # Si factor es "reducción", el factor de zoom es 1/reducción.
    # Mantengamos factor como el cambio de tamaño físico (factor > 1 agranda)
    zoom_factor = 1.0 / factor 

    resampled_data = zoom(input_data, zoom_factor, order=order)
    return resampled_data

# --- Código de Prueba (Temporal) ---
# Esto nos permite probar el módulo directamente, 
# reemplazando el script antiguo.
if __name__ == '__main__':
    imagen_fits = '../../pruebas.fits' # Asume que se ejecuta desde src/cmod/
    data = fits.getdata(imagen_fits)

    resize_factor = 1.5 # Factor para hacerla 1.5x más pequeña

    data_reducida = geometric_resample(data, resize_factor, order=1)

    fits.writeto('../../pruebas_resize_new.fits', data_reducida, overwrite=True)
    print("Imagen resampleada guardada como 'pruebas_resize_new.fits'")