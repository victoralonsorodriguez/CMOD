from astropy.io import fits
import numpy as np
from scipy.ndimage import zoom

# Cargar imagen FITS
imagen_fits = 'pruebas.fits'  # Cambia por tu archivo
data = fits.getdata(imagen_fits)

# Factor de reducción (ejemplo: reducir a 1.5 veces más pequeña)
factor = 1.5
factor_resize = 1 / factor  # Para hacer la imagen más pequeña

# Reducir imagen con interpolación
data_reducida = zoom(data, factor_resize, order=1)  # 'order=1' usa interpolación bilineal

# Guardar la nueva imagen FITS
fits.writeto('pruebas_resize.fits', data_reducida, overwrite=True)

print(f"Imagen reducida guardada como 'imagen_reescalada.fits'")