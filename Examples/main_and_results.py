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
   
    
#RESULTS OF PRINT
 
# moshits_Q2_atm01_0014_10PeV_15_000_c010
# done
# duration =  550.0 ns
# axis x, (mm) y (mm), d (m)): (-72.29147311947149, 109.50821693614641, 231.07075847659786)
# theta, phi, a0, a1, a2: (0.26836619222268515, 4.919588940193713, 85.50520032459649, -0.32077714631746973, 0.004931077593639349)
     
# moshits_Q2_atm01_0014_10PeV_15_000_c001
# done
# duration =  525.0 ns
# axis x, (mm) y (mm), d (m)): (-97.51183246086578, -17.66539678704814, 179.2821962874435)
# theta, phi, a0, a1, a2: (0.25144340994796516, 4.955215626795071, 75.97296984135814, -0.19404964710314299, 0.004412811954492142)
     
# moshits_Q2_atm01_0014_10PeV_15_000_c006
# done
# duration =  750.0 ns
# axis x, (mm) y (mm), d (m)): (3.689519982780385, -151.032113050369, 261.32259903635475)
# theta, phi, a0, a1, a2: (0.21984899401933483, 4.952921377402563, 72.67528334551473, -0.3006751370177325, 0.004878761181191166)
     
# moshits_Q2_atm01_0014_10PeV_15_001_c007
# done
# duration =  625.0 ns
# axis x, (mm) y (mm), d (m)): (132.10245124346687, 29.82750324820622, 241.50695643250214)
# theta, phi, a0, a1, a2: (0.27241804378980045, 3.0845345337925947, 81.7491014014713, -0.20596210136696275, 0.004383255463439998)
     






    