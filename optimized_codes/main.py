%%time

### =========================================================================
### =========================================================================

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


### =========================================================================
### =========================================================================

## === !!! CHANGE ONLY HERE !!!

telescope='SPHERE-2'        #Change: name of the telescope  - SPHERE-2 / SPHERE-3
type_analyse='electronic'   #Change: type of analyse - only_sig / electronic / experiment

if type_analyse!='experiment': 
    En, nuclei, H = 10, 'Fe', 900 #Change: energy, nuclei, altitude
else:
    En, nuclei, H = '', '', ''

maindir = f'/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/{telescope}/{type_analyse}/'  #Change the directory path
if type_analyse!='experiment':
    maindir=maindir+f'{type_analyse}_{En}PeV_{nuclei}_{H}m/' #Change the name of the sample

## === !!! END CHANGE ONLY HERE !!!

## =========================================================================
### =========================================================================

## === DASHBOARD

for filename in os.listdir(maindir):
    if type_analyse=='experiment':
        if not filename.endswith('.txt'):
            continue
    #Initial values
    initial_values=Initial_values(type_analyse, telescope, maindir+filename,nuclei,En,H)
    #Open files
    files=Openfiles(initial_values)
    #Algorithm of angles determination
    algorithm=Functions_algorithm(initial_values,files)
    algorithm.length
    
    # #Condition for events selection
    conditions=Conditions(algorithm,initial_values)
    
    #Results and Figures
    results_pictures=Results(initial_values,algorithm)
    # results_pictures.fig_sum_impulse #Fig
    # results_pictures.front  #Fig
    results_pictures.angles #Results
   
    

    