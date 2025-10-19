#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script principal para ejecutar el pipeline de análisis CMOD.

Este script ejecuta las dos fases principales:
1. Análisis Científico: Procesa los FITS y genera los CSVs.
2. Ploteo: Lee los CSVs y genera los plots finales del paper.
"""

import sys
from cmod.utils import Cronometro
from cmod.pipeline import run_analysis_pipeline
from cmod.plotting import run_plotting_pipeline


if __name__ == '__main__':

    # (Opcional: aquí podríamos tener un argparse simple 
    # para decidir si correr 'analysis', 'plot', o 'all')

    print("--- Iniciando Pipeline CMOD ---")

    cronometro_total = Cronometro()
    cronometro_total.iniciar()

    #######################################################
    ###------------------CMOD ANALYSIS------------------###
    #######################################################

    print("\n[FASE 1/2] Ejecutando el pipeline de ANÁLISIS...")
    try:
        run_analysis_pipeline()
        print("[FASE 1/2] Pipeline de ANÁLISIS completado.")
    except Exception as e:
        print(f"¡Error fatal en el pipeline de ANÁLISIS!:\n{e}")
        sys.exit(1) # Salir si el análisis falla


    ##########################################################
    ###------------------PLOTTING RESULTS------------------###
    ##########################################################

    print("\n[FASE 2/2] Ejecutando el pipeline de PLOTEO...")
    try:
        run_plotting_pipeline()
        print("[FASE 2/2] Pipeline de PLOTEO completado.")
    except Exception as e:
        print(f"Error en el pipeline de PLOTEO:\n{e}")

    print("\n--- Pipeline CMOD Finalizado ---")
    cronometro_total.detener()