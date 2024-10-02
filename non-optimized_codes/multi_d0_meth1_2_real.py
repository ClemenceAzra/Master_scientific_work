%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import glob
import os
import random

###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

#Type of telescope
name_telesc='SPHERE-2' #SPHERE-2, SPHERE-3

#Altitude of the telescope
H=500 #500, 900, 1000 m

#Type of analyse
type_analyse='only_sig'
# type_analyse='sig_&_bg'
# type_analyse='electronic'

#Energy
En=30

###=========================================================================

##=== FONCTIONS

#===============

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{}PeV_P_{}m'.format(En,H) #Signal only

    def read_mosaic_hits_file(path, step = 12.5):
        
        a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+')     
        a[4]=a[4]-min(a[4]) # номер столбца со временем - 4
       
        pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[4])}) #creation of DataFrame to facilite calculations
        q = np.zeros((109, 1020))  # empty array109 PMT * 1020 cells
        cells_here=np.arange(0,1021*step,step)
     
        shift = random.randint(400, 600) #first bin of the signal
        
        for i in num_pmt:
            impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
            n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
            n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
            q[i]=n_phot_shift

#Signal with background

            if type_analyse=='sig_&_bg':
                data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

                q[i]=n_phot_shift+np.random.choice(data_photons_bg['{}'.format(i)].dropna(),size=1020)
                
        return pd.DataFrame(q.T)
#Electronics

if type_analyse=='electronic':
    
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_P_{}m'.format(En,H) #Electronic signal
    
    def read_electonic_event_file(path):
        event = pd.read_csv(path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip')
        for pmt in range(event.shape[1]):
            event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
            event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
        return event


name_files = glob.glob(os.path.join(path,'*'))

#== Start values for multistart

def theta_initt(theta_init):
    return theta_init*len(theta_init)

#== Fonction of the approximation of the front

def tau_fnc(theta,phi,a0,a1,a2):
    x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
    y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
    z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
    R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
    return a0+a1*R+a2*R**2+z_casc/c_ns

###=========================================================================

##Coordinates of PMT in mosaic

#PMT on mosaic: x (mm),y (mm), num
if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/mosaic_pmt_coords_df.csv')

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel

###=========================================================================

##=== CONSTANT VALUES

cells=np.arange(0,1020*12.5,12.5) #1020 cells with 12.5 ns each
num_pmt=np.arange(0,109) #Num of PMT from 1 to 109 

if name_telesc=='SPHERE-2':
    if H==500:
        a=0.0005444285453203233 #a, b for different H --> find in another work
        b=0.9051146902370498
    elif H==900:
        a=0.0008982897531551649
        b=1.6455812746636922
       
if name_telesc=='SPHERE-3':
    b=(1.1314175328971585)/(1000/H) #If H=500 - a --> a/2


c_ns=0.3 #speed of light m/ns

part_of_int=0.5 #part of integral 

#=====FNC

#a0 fnc
# a_a0=-0.630484 #more in the article "Алгоритм восстановления направления прихода ШАЛ для телескопа СФЕРА"
# b_a0=53.541113661

#Multistart
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

# real_theta_deg=0.1745*180/np.pi
# real_phi_deg=3.0564*180/np.pi

theta_real=0.3277
phi_real=3.0564

# # theta_real=0.2175
# # phi_real=2.982

# # theta_real=0.1745
# # phi_real=3.0564

##=== PRELIMINARY CALCULATIONS

#==Translation x,y from mos to snow

if name_telesc=='SPHERE-2':
    x_snow=-b*x_y_mos['x']
    y_snow=-b*x_y_mos['y']

if name_telesc=='SPHERE-3':
    x_snow=-b*np.array(pixel_data.groupby('segment').mean()['x'])
    y_snow=-b*np.array(pixel_data.groupby('segment').mean()['y'])

#==Translation t from mos to snow (path from mos to snow)
  
t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns

###=========================================================================

real_x0=[]
real_y0=[]

max_x0=[]
max_y0=[]

grav_center_x0=[]
grav_center_y0=[]
#===
for files in range(0,len(name_files)):
    print(name_files[files])
    #Open file
    if type_analyse!='electronic':
        event = read_mosaic_hits_file(name_files[files])
    if type_analyse=='electronic':
        event = read_electonic_event_file(name_files[files])
    
    if len(event)<1020:
        continue
    
    ###========================================================================
    
    
    #Search the impulse in the отклик
    
    pulse = event.sum(axis=1) # total impulse in 1020 cells
    
    if max(pulse[:-1])<max(pulse[0:400])*2:
        continue
    
    imax = pulse[:-1].idxmax() #index of the max impulse --> delete the last cell -> > max of impulse
    start, stop = imax, imax
    
    maxpart = 1
    bgd = maxpart * max(pulse[0:400])
    while pulse[start] > bgd:
        start -= 1
    while pulse[stop] > bgd:
        stop += 1
    
    margin_left,margin_right= 30, 30 
    
    stop  += margin_right
    start -= margin_left
    
    event_diapason=event[start:stop] #diapason of the impulse - n photons
    cells_diapason=cells[start:stop] #diapason of the impulse - time
    
    ###============================================
    
    #Save the PMT and bins in with no BG 
    
    saved_pmt=[]
    t_front_reconstr=[]
    x_front_reconstr=[]
    y_front_reconstr=[]
    number_photons_i=[]
    
    impulse_each_pmt=[]
    cells_each_pmt=[]
    
    for i in range(0,event_diapason.shape[1]):
        event_pmt=event_diapason[i]
        
        #Condition 1 : delete the PMT with max signal < max BG
        if max(event_pmt)<max(event[i][0:400]):
            continue
        
        imax_pmt = event_pmt.idxmax() #index of the max impulse
        start_pmt, stop_pmt = imax_pmt, imax_pmt

        bgd_pmt = max(event[i][0:400])
        while event_pmt[start_pmt] > bgd_pmt and start_pmt>event_pmt.index[0]:
            start_pmt -= 1  
            
        while event_pmt[stop_pmt] > bgd_pmt and stop_pmt<event_pmt.index[-1]:
            stop_pmt += 1
            
        # #Impulse in the finded diapason
        impulse=event_pmt[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
        cells_impulse=cells_diapason[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]

        #Condition 2: delete bins < max BG
        impulse_filtrate=impulse[impulse>max(event[i][0:400])] 
        cells_impulse_filtrate=cells_impulse[impulse>max(event[i][0:400])] 

            
        #Condition 3 : delete the PMT with 0 bins
        if len(cells_impulse_filtrate)<1: #if number of bins in PMT less than 4 - low signal
            continue
        
        #Aoment of the arrival of shower in the PMT
        t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

        #SAVE DATA
        
        t_front_reconstr.append(t_arrival-t_path[i]) #t - obtain front
        x_front_reconstr.append(x_snow[i]) #x - obtain front
        y_front_reconstr.append(y_snow[i]) #y - obtain front
        
        saved_pmt.append(num_pmt[i]) # Nº of the saved PMT with signal
        number_photons_i.append(impulse_filtrate.sum()) #number of signal photons in each PMT

    ###============================================
    
    ###============================================
    
    if len(x_front_reconstr)<4:
        continue
    
    #Axis of the shower

    x0_front=x_front_reconstr[number_photons_i.index(max(number_photons_i))] 
    y0_front=y_front_reconstr[number_photons_i.index(max(number_photons_i))]
    t0_front=t_front_reconstr[number_photons_i.index(max(number_photons_i))]
    
    d0_from_center=np.sqrt(x0_front**2+y0_front**2)

    d_x_to_x_max=np.sqrt( (np.array(x_front_reconstr)-x0_front)**2  + (np.array(y_front_reconstr)-y0_front)**2  )

    if d0_from_center>180:
            continue
        
    elif d0_from_center>140 and d0_from_center<180:
        x_circle_around_max=np.array(x_front_reconstr)[np.where(d_x_to_x_max<50)[0]]
        y_circle_around_max=np.array(y_front_reconstr)[np.where(d_x_to_x_max<50)[0]]
        number_photons_around_max=np.array(number_photons_i)[np.where(d_x_to_x_max<50)[0]]

    elif d0_from_center<140:
        x_circle_around_max=np.array(x_front_reconstr)[np.where(d_x_to_x_max<100)[0]] 
        y_circle_around_max=np.array(y_front_reconstr)[np.where(d_x_to_x_max<100)[0]]
        number_photons_around_max=np.array(number_photons_i)[np.where(d_x_to_x_max<100)[0]]

    x0_grav_center=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
    y0_grav_center=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)

    x0_real=float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[2])
    y0_real=float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[3])
    
    #SAVE DATA

    real_x0.append(x0_real)
    real_y0.append(y0_real)
    
    max_x0.append(x0_front)
    max_y0.append(y0_front)

    grav_center_x0.append(x0_grav_center)
    grav_center_y0.append(y0_grav_center)

res_d=np.array([real_x0,real_y0,max_x0,max_y0,grav_center_x0,grav_center_y0])
