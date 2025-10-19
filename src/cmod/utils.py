import numpy as np
import os
import threading
import shutil
import atexit
import time
import sys

class Cronometro:
    def __init__(self):
        """Inicializa el cronómetro pero no lo inicia."""
        self.inicio = None
        self.cronometro_activo = False
        atexit.register(self.detener)  # Se detendrá automáticamente al final

    def iniciar(self):
        """Inicia el cronómetro en un hilo separado."""
        if self.cronometro_activo:
            return  # Evita iniciar múltiples veces

        self.inicio = time.time()
        self.cronometro_activo = True
        self.hilo = threading.Thread(target=self._ejecutar, daemon=True)
        self.hilo.start()

    def _ejecutar(self):
        """Muestra el cronómetro en la última línea de la terminal mientras está activo."""
        while self.cronometro_activo:
            tiempo_transcurrido = time.time() - self.inicio
            tiempo_formateado = self._formatear_tiempo(tiempo_transcurrido)
            columnas = shutil.get_terminal_size().columns  # Obtener el ancho de la terminal
            
            # Mover el cursor a la última línea y sobrescribir
            sys.stdout.write(f"\033[{shutil.get_terminal_size().lines}E")  
            sys.stdout.write(f"\033[KTime running: {tiempo_formateado}".ljust(columnas) + "\r")
            sys.stdout.flush()
            time.sleep(1)

    def detener(self):
        """Detiene el cronómetro y muestra el tiempo total."""
        if not self.cronometro_activo:
            return

        self.cronometro_activo = False
        tiempo_final = time.time() - self.inicio
        tiempo_formateado = self._formatear_tiempo(tiempo_final)
        columnas = shutil.get_terminal_size().columns

        # Mostrar el tiempo total al finalizar
        sys.stdout.write(f"\033[{shutil.get_terminal_size().lines}E")  
        sys.stdout.write(f"\033[KTotal time: {tiempo_formateado}\n".ljust(columnas))
        sys.stdout.flush()

    def _formatear_tiempo(self, segundos):
        """Convierte el tiempo en segundos a formato hh:mm:ss"""
        horas = int(segundos // 3600)
        minutos = int((segundos % 3600) // 60)
        segundos = int(segundos % 60)
        return f"{horas:02}:{minutos:02}:{segundos:02}"

def create_folder(folder_path,
                  overwrite = False):
    
    # Check if the folder already exists
    if os.path.isdir(folder_path) == True:
        
        # If overwrite is selected
        if overwrite == True:
            shutil.rmtree(folder_path)
            os.mkdir(folder_path)

        # If not, pass
        else:
            pass
    
    # If the folder is not already created
    else:
        os.mkdir(folder_path)

def round_number(num,dec):
    
    if isinstance(num, (int, float, np.float32, np.float64)):
        num = round(num,dec)
        
    elif isinstance(num, np.ndarray):
        num = np.round(num,dec)   

    return num


###------------------ELLIPTICITY AND AXIS RATIO------------------###
####################################################################

def ell_to_axrat(ell):
    axrat = 1 - ell
    return round_number(axrat,3)

def axrat_to_ell(axrat):
    ell = 1 - axrat
    return round_number(ell,3)

###------------------RADIANS AND DEGREES------------------###
#############################################################

# Changing from degrees to radian
def deg_to_rad(deg):
    
    rad = (deg) * (np.pi / 180)

    return round_number(rad,3)

def rad_to_deg(rad):
    
    deg = (rad) * (180 / np.pi)
        
    return round_number(deg,3)


def rad_to_deg_abs(rad):
    
    deg = ((rad) * (180 / np.pi)%360)
        
    return round_number(deg,3)