### IMPORT MODULES

import pandas as pd
import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt
import os
from initial_values import *
from open_files import *
from functions_algorithm import *
from conditions import *
from results import *
import warnings
warnings.filterwarnings('ignore')

### =========================================================================
### =========================================================================

## === !!! CHANGE ONLY HERE !!!

telescope='SPHERE-3'        #Name of the telescope  -> SPHERE-2 / SPHERE-3
type_analyse='only_sig'   #Type of analyse -> only_sig / electronic / experiment

En = 10 #Energy, PeV
nuclei = 'P' #N, Fe
H = 500 #detector altitude, m

maindir = f'/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/{telescope}/{type_analyse}/{type_analyse}_{En}PeV_{nuclei}_{H}m/'  #Directory of the event files
dir_params = '/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/Initial_data/' #Directory of initial parameters

## === !!! END CHANGE ONLY HERE !!!

## =========================================================================
### =========================================================================

## === DASHBOARD

for filename in os.listdir(maindir):
   
    #Do not read the file "DS store"
    
    if type_analyse == 'experiment':
        if not filename.endswith('.txt'):
            continue
    if type_analyse != 'experiment':
        if not filename[-1].isdigit(): 
            continue
    
    #Pre opening, obligatory!
    
    initial_values=Initial_values(telescope, type_analyse, En, nuclei, H, maindir+filename, dir_params)
    files=Openfiles(initial_values)
    fnc_algorithm=Functions_algorithm(initial_values, files)
    
    conditions=Conditions(fnc_algorithm,initial_values) 
    if conditions.algorithm_condition != 'done':
        continue
    
    #Results of the functions, choose which result see
    
    print(filename)
    print('duration = ', fnc_algorithm.length(), 'ns')
    print('axis x, (mm) y (mm), d (m)):',fnc_algorithm.axis)
    print('theta (rad), phi (rad), a0, a1, a2:', fnc_algorithm.angles)
    print('     ')
    
    #Results and Figures
    results_pictures=Results(initial_values, fnc_algorithm)
    # results_pictures.fig_sum_impulse #Fig
    # results_pictures.front  #Fig
    # print(results_pictures.angles) #Results
   
    

    