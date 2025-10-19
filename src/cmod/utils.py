import numpy as np
import os
import re
import threading
import shutil
import atexit
import time
import sys
from datetime import datetime



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




def version_directory(cwd,base_name,original_dir='.',create_dir=False):

    pattern = re.compile(rf'{re.escape(base_name)}_V(\d+)')
    versions = []

    # Searching for version
    for dir_name in os.listdir(cwd):
        match = re.match(pattern, dir_name)
        if match:
            versions.append(int(match.group(1)))
    
    # If no version is found
    if len(versions)==0:
        versions.append(0)
    
    # Determining the new version
    versions = np.array(versions)
    previus_version = max(np.nanmax(versions),0)
    next_version = previus_version + 1
    
    # Creating the folder
    next_dir_path = f'{original_dir}_V{next_version}'
    if create_dir == True:
        create_folder(next_dir_path)
    
    return int(previus_version)


def version_file(galaxy_name, version, analysis_version, psf_size,archivo,analysis_programe):

    datetime_current = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    date, time = datetime_current.split(' ')
    
    nueva_linea = (f'VERSION={version} PROGRAMME={analysis_programe} DATE={date} TIME={time} ' 
                   f'ANALYSIS_MODE={analysis_version} PSF={psf_size}\n')

    # Comprobar si el archivo ya existe
    existe = os.path.exists(archivo)

    with open(archivo, 'a' if existe else 'w') as f:
        # Si el archivo no existe, escribir el encabezado
        if not existe:
            f.write(f'#-------GALAXY {galaxy_name}-------#\n\n')
        f.write(nueva_linea)


def version_file_last(file):

    if not os.path.exists(file):   
        return 0

    versions = []
    with open(file, 'r') as f:
        for line in f:
            match = re.search(r'VERSION=(\d+)', line)
            if match:
                versions.append(int(match.group(1)))
    
    max_version = max(versions)

    return int(max_version)