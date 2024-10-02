%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import glob
import os


###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

#Type of telescope
name_telesc='SPHERE-2' #SPHERE-2, SPHERE-3

#Altitude of the telescope
H=500 #500, 900, 1000 m

#Type of analyse
type_analyse='only_sig'
# type_analyse='sig_&_bg'
type_analyse='electronic'
# type_analyse='sph3'
# type_analyse='experiment'

#Energy
En=30

part_of_int=0.5

nuclei='pro'

###=========================================================================
###=========================================================================

##=== OPEN FILES

#===============
#== For all types of analyse

#PMT on mosaic: x (mm),y (mm), num
if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/mosaic_pmt_coords_df.csv')
    total_number_of_pmt=109
    #Open the number of the PMT around a PMT
    circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/circle_of_pmt_around_pmt.csv',)#===============
    coeff_amp=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/coeff_amp.csv',header=None,sep='\s+')

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel
    total_number_of_pmt=379
    
    circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/circle_of_pmt_around_pmt_sph3.csv',)#===============

# #===============

# #Open the number of the PMT around a PMT

# circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/circle_of_pmt_around_pmt.csv',)#===============

#===============

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{}PeV_P_{}m'.format(En,H) #Signal only
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res_500m'.format(En,H) #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_{}_{}m'.format(En,nuclei,H) #Electronic signal #044 #072
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/little'.format(En,nuclei,H) #Electronic signal #044 #072

    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res_{}m'.format(H) #Electronic signal d0=147 m
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res/mosaic_hits_m01_Fe_10PeV_10-20_038_c005' #Electronic signal #044 #072

#Experiment

if type_analyse=='experiment':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/3_experiment'

#======

#SPHERE 3

if name_telesc=='SPHERE-3':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/sph3/{}m_{}PeV_{}'.format(H,En,nuclei)
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/sph3/{}m_{}PeV_{}'.format(H,En,nuclei)

#Open files for real angles


if name_telesc=='SPHERE-2':
    real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/angles_theta',header=None,sep='\s+')[0]) #number of PMT
    real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/4_data_init/angles_phi',header=None,sep='\s+')[0]) #number of PMT
    theta_real=real_theta[0]
    phi_real=real_phi[0]

elif name_telesc=='SPHERE-3':
    real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_{}'.format(nuclei),header=None,sep='\s+')[0]) #number of PMT
    theta_real=15/180*np.pi
    phi_real=real_phi[1]
    



name_files = glob.glob(os.path.join(path,'*'))

###=========================================================================
###=========================================================================

##=== FUNCTION

##===================

#=== READ FILES

# Return event with the type of analyse
if name_telesc=='SPHERE-2':
    num_t_event=4
    lim_search=100
if name_telesc=='SPHERE-3': 
    num_t_event=5
    lim_search=50


# Only signal
def read_mosaic_hits_file(path, step = 12.5):

    a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+') 
    a[num_t_event]=a[num_t_event]-min(a[num_t_event])  
    pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[num_t_event])}) #creation of DataFrame to facilite calculations
    q = np.zeros((total_number_of_pmt, 1020))  # empty array109 PMT * 1020 cells
    cells_here=np.arange(0,1021*step,step)
       
    shift = 500 #first bin of the signal
    
    for i in num_pmt:
        impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
        n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
        n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
        q[i]=n_phot_shift

    return pd.DataFrame(q.T)

# ===

#Signal + background
def read_mosaic_hits_file_background(path, step = 12.5):

    a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+')     
    a[num_t_event]=a[num_t_event]-min(a[num_t_event]) # номер столбца со временем - 4
       
    pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[num_t_event])}) #creation of DataFrame to facilite calculations
    q = np.zeros((total_number_of_pmt, 1020))  # empty array109 PMT * 1020 cells
    cells_here=np.arange(0,1021*step,step)
       
    shift = 500 #first bin of the signal
    
    for i in num_pmt:
        impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
        n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
        n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
        q[i]=n_phot_shift    
        q[i]=n_phot_shift+np.random.choice(data_photons_bg['{}'.format(i)].dropna(),size=1020)
            
    return pd.DataFrame(q.T)

# ===

#Electronic for SPHERE-2
def read_electonic_event_file(path):
    event = pd.read_csv(path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip')
    for pmt in range(event.shape[1]):
        event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
        event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
    for i in range(event.shape[1]):
        event[i]=event[i]*coeff_amp[1].iloc[i]
    return event


#Electronic for SPHERE-3
def read_electonic_event_file_sph3(path):
    event = pd.read_csv(path, header=None, sep='\s+',on_bad_lines='skip').T
    event['pixel'],event['segment']=event.index%8,event.index//8 
    event=event[event['pixel']!=7].groupby('segment').mean().drop(['pixel'], axis=1).T  #delete each 8th pixel and average each segment + #drop column pixel 
    for pmt in range(event.shape[1]):
        event[pmt] -= event[0:100][pmt].mean() #substraction of the piedestal
    return event


# ===

#Experiment for SPHERE-2

def read_experiment_sph2(path):
    event = pd.read_csv(path,skiprows=40,skipfooter=2144,engine='python',sep='\s+',header=None)
    event=event.drop(columns=[0,110,111,112])
    event.columns=event.columns-1
    for pmt in range(event.shape[1]):
        event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
        event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
    for i in range(event.shape[1]):
        event[i]=event[i]*coeff_amp[1].iloc[i]
    return event


# _____________________________


# Return start values of multistart

def theta_initt(theta_init):
    return theta_init*len(theta_init)

##===================

#=== DIAPASON OF THE IMPULSE

def diapason_impulse(event):
    pulse = event.sum(axis=1) # total impulse in 1020 cells
    imax = pulse[:-1].idxmax() #index of the max impulse --> delete the last cell -> > max of impulse
    start, stop = imax, imax
    
    maxpart = 1
    background=np.mean(pulse[1:pied])+3*np.std(pulse[1:pied])
    bgd = maxpart * background
    
    while pulse[start] > bgd and start!=pulse.index[0]:
        start -= 1
    while pulse[stop] > bgd and stop!=pulse.index[-1]:
        stop += 1
    
    margin_left,margin_right= 30, 30 
    
    stop  += margin_right
    start -= margin_left
    
    event_diapason=event[start:stop] #diapason of the impulse - n photons
    cells_diapason=cells[start:stop] #diapason of the impulse - time
    
    return event_diapason,cells_diapason

##===================

#=== AXIS OF THE SHOWER

def axis(x0_max,y0_max,all_window):
    
    d_pmt_to_d_max=np.sqrt( (np.array(x_mos)-x0_max)**2  + (np.array(y_mos)-y0_max)**2  ) #distance of each PMT from the max PMT, m
    
    if d0_from_center>d0_before_c and d0_from_center<d0_center: #if the max PMT is on the circle before the last one
        max_d=d_six #the minimum distance of the neighors around
    
    elif d0_from_center<d0_before_c: #if the max PMT is inside
        max_d=d_double_six #the minimum distance of the neighors around
        
    number_photons=list(all_window.sum()) #max of each PMT
    x_circle_around_max=x_mos[np.where(d_pmt_to_d_max<max_d)[0]] #x around the max PMT
    y_circle_around_max=y_mos[np.where(d_pmt_to_d_max<max_d)[0]] #y around the max PMT
    number_photons_around_max=np.array(number_photons)[np.where(d_pmt_to_d_max<max_d)[0]] #max PMT around the max PMT
    
    x0=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max) #x by gravity center
    y0=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max) #y by gravity center
    return x0, y0

##===================

#=== TIME OF ARRIVAL

#Return PMT than can be analysed

def saved_PMT(df_front_N_1,df_front_N_2):
    for k in range(0,20):
        good_pmt=[]
        bad_pmt=[]
        pmt_without_ngbr=[]
        for i in range(0,len(df_front_N_2)):
            looked=df_front_N_2.iloc[i] #number of the PMT in the table 2
            PMT_around=circle_of_pmt_around_pmt.loc[looked.name].dropna()[1:]  #PMT around the PMT, delete the central PMT
            
            PMT_around_in_table_1=PMT_around[PMT_around.isin(df_front_N_1['PMT'])]
            PMT_around_in_table_2=PMT_around[PMT_around.isin(df_front_N_2['PMT'])]
             
            if len(PMT_around_in_table_1)==0 and len(PMT_around_in_table_2)==0:
                pmt_without_ngbr.append(looked['PMT'])
                continue
            
            if len(PMT_around_in_table_1)==0: #No neighbours of the PMT of table 2 in the table 1
                continue
            
            else:
                mean_time=df_front_N_1.loc[PMT_around_in_table_1]['t_max'].mean() #mean of the sure PMT around the examinated PMT
                
                if looked['t_max'] <= mean_time+lim_search and  looked['t_max'] >= mean_time-lim_search:
                    good_pmt.append(looked['PMT'])
                else:
                    bad_pmt.append(looked['PMT'])
                    continue
        
        df_front_N_1=pd.concat([df_front_N_1,df_front_N_2.loc[good_pmt]])
        df_front_N_2=df_front_N_2.drop(good_pmt+bad_pmt+pmt_without_ngbr) #Delete sure PMT
     
    return df_front_N_1

#Delete PMT if there is less than 3 PMT around
def delete_less_3_nbr(good_PMT):
    bad_PMT=[]
    for i in range(0,len(good_PMT)):
        if len(good_PMT['PMT'][good_PMT['PMT'].isin(circle_of_pmt_around_pmt.loc[good_PMT.iloc[i][0]].dropna())].drop(good_PMT.iloc[i][0]))<2:
            bad_PMT.append(int(good_PMT.iloc[i][0]))
    
    good_PMT_new=good_PMT.drop(bad_PMT)

    return good_PMT_new

#Return time of arrival and N photons by part of integral

def t_part_of_int(PMT_total,index_total):
    t_front_reconstr=[]
    N_front_reconstr=[]
    x_front_reconstr=[]
    y_front_reconstr=[]
    PMT=[]
    for i in range(0,len(x_front)):
        # print(i)
        window=all_window[PMT_total[i]]
        event_pmt=window
    
        imax_pmt = int(index_total[i]) #index of the max impulse in the diapason
        start_pmt, stop_pmt = imax_pmt, imax_pmt
         
        if type_analyse=='only_sig':
            threshold=0
        else:
            threshold=mean_std[PMT_total[i]]

        
        while event_pmt[start_pmt] > threshold and start_pmt>event_pmt.index[0]:
            start_pmt -= 1  
             
        while event_pmt[stop_pmt] > threshold and stop_pmt<event_pmt.index[-1]:
            stop_pmt += 1
             
        # #Impulse in the finded diapason
        # impulse_filtrate=event_pmt[start_pmt+1-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
        # cells_impulse_filtrate=cells_diapason[start_pmt+1-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
    
        # impulse_filtrate=event_pmt[start_pmt+1-min(event_pmt.index):imax_pmt+1]
        # cells_impulse_filtrate=cells_diapason[start_pmt+1-min(event_pmt.index):imax_pmt+1]

        impulse_filtrate=event_pmt[imax_pmt+1:stop_pmt-min(event_pmt.index)]
        cells_impulse_filtrate=cells_diapason[imax_pmt+1:stop_pmt-min(event_pmt.index)]

        if len(impulse_filtrate)<1:
            continue
    #___________________
    
    #Aoment of the arrival of shower in the PMT
        # t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]
    
        t_arrival=cells_impulse_filtrate[len(impulse_filtrate[impulse_filtrate<impulse_filtrate.max()/2])-1]

        t_front_reconstr.append(t_arrival-t_path[int(PMT_total[i])])
        N_front_reconstr.append(impulse_filtrate.sum())
        x_front_reconstr.append(x_snow[PMT_total[i]])
        y_front_reconstr.append(y_snow[PMT_total[i]])
        PMT.append(PMT_total[i])

    return t_front_reconstr, N_front_reconstr, x_front_reconstr, y_front_reconstr, PMT

##===================

#=== ANGLES

def angles(x_front,y_front,t_fnc,N_photons):
    
    def fnc(theta,phi,a0,a1,a2):#определение функции с не заданными параметрами
            x_casc=((np.cos(theta)*np.cos(phi))*(x_front)+(np.cos(theta)*np.sin(phi))*(y_front)) #координаты x фотонов в системе ливня
            y_casc=(-np.sin(phi)*(x_front)+np.cos(phi)*(y_front)) #координаты y фотонов в системе ливня
            z_casc=((np.sin(theta)*np.cos(phi))*(x_front)+(np.sin(theta)*np.sin(phi))*(y_front)) #координаты z фотонов в системе ливня
            R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
            tau=a0+a1*R+a2*R**2+z_casc/c_ns #аппроксимированое время
            s=(((t_fnc-tau)**2))*N_photons #функцию которую надо аппроксимировать
            Smin=np.sum(s)
            return Smin
        
    theta_multi_start=[]
    phi_multi_start=[]
    min_s=[]
        
    a0_multi_start=[]
    a1_multi_start=[]
    a2_multi_start=[]
    
    for th in range(0,len(theta_initt)):
            param_fix=Minuit(fnc, theta=theta_initt[th],phi=phi_initt[th],a0=1,a1=1,a2=1)
            param_fix.limits['theta']=(0,60/180*np.pi)
            param_fix.limits['phi']=(0,2*np.pi)
    
            param_fix.migrad()
            theta_multi_start.append(param_fix.values[0]) 
            phi_multi_start.append(param_fix.values[1])
            
            a0_multi_start.append(param_fix.values[2]) 
            a1_multi_start.append(param_fix.values[3])
            a2_multi_start.append(param_fix.values[4])
            
            min_s.append(param_fix.fval)
        
    theta=np.array(theta_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    phi=np.array(phi_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    
    a0=np.array(a0_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    a1=np.array(a1_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    a2=np.array(a2_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    
    return theta, phi, a0, a1, a2, param_fix


###=========================================================================
###=========================================================================

##=== CONSTANT VALUES


# === CELLS AND PIEDESTAL

cells=np.arange(0,1020*12.5,12.5) #Length of the electronic response, bin=12.5 ns *1200
pied=400 #piedestal
num_pmt=np.arange(0,total_number_of_pmt) #number of PMT for 0 to the last

if type_analyse=='electronic' and name_telesc=='SPHERE-3': #different cells and piedestal
    pied,cells=100,np.arange(0,500*12.5,12.5)

if type_analyse=='experiment' and name_telesc=='SPHERE-2': #different cells and piedestal
    pied,cells=200,np.arange(0,915*12.5,12.5)
    
    
# === Parameters

if name_telesc=='SPHERE-2':
    
    #Coeff for translate x mos --> x snow
    b=0.9051146902370498*(H/500) 
    
    #Values for the axis
    d0_center=192 #last circle on mosaic
    d0_before_c=154 #before last circle on mosaic
    d_six=50  #distance PMT and 1 circle around
    d_double_six=100 #distance PMT and 2 circles around
    
    #Diapason for search the impulse by BFS
    lim_search=100 #ns


       
if name_telesc=='SPHERE-3':
    
    #Coeff for translate x mos --> x snow
    b=0.5657087664485793*(H/500) #Coeff for translate x mos --> x snow
     
    #Values for the axis
    d0_center=292  #last circle on mosaic
    d0_before_c=262 #before last circle on mosaic
    d_six=32 #distance PMT and 1 circle around
    d_double_six=70  #distance PMT and 2 circles around
    
    #Diapason for search the impulse by BFS
    lim_search=50 #ns

# === PHYSICAL VALUE

c_ns=0.3 #speed of light m/ns


#=====FNC

# a2=0.0005

#Paramatrization
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

real_theta_deg=0.1745*180/np.pi
real_phi_deg=3.0564*180/np.pi

# theta_real=0.3277 #[0]
# phi_real=3.0564 #[0]

# theta_real= 0.2718 #[files 008]
# phi_real=2.4855 #[files 008]

# theta_real=0.3332 #[14]
# phi_real=5.7728 #[14]

# theta_real=0.2175
# phi_real=2.982

# theta_real=0.1745
# phi_real=3.0564

##=== PRELIMINARY CALCULATIONS

#==Translation x,y from mos to snow

if name_telesc=='SPHERE-2':
    x_mos=np.array(x_y_mos['x'])
    y_mos=np.array(x_y_mos['y'])
    x_snow=-b*x_mos
    y_snow=-b*y_mos

if name_telesc=='SPHERE-3':
    x_mos=np.array(pixel_data.groupby('segment').mean()['x'])
    y_mos=np.array(pixel_data.groupby('segment').mean()['y'])
    x_snow=-b*x_mos
    y_snow=-b*y_mos


#==Translation t from mos to snow (path from mos to snow)
  
t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns

###========================================================================
###========================================================================

delta_by_max=[]
delta_by_part_of_int=[]
files_saved=[]
total_nbr_pmt=[]
lenght_impulse_total=[]
d0_all=[]
big_max_summed_sig=[]
event_diapason_tot=[]
theta_max=[]
phi_max=[]
x0_real=[]
y0_real=[]

x0_max_tot=[]
y0_max_tot=[]

x0_grav=[]
y0_grav=[]

df_results_all=[]

for files in range(0,len(name_files)):
    # print(name_files[files])

    #Open file
    # Only signal
    if type_analyse=='only_sig':
        event = read_mosaic_hits_file(name_files[files])
        
    # Signal with the background
    if type_analyse=='sig_&_bg':
        event = read_mosaic_hits_file_background(name_files[files])

    #Signal with the electronic
    if type_analyse=='electronic' and name_telesc=='SPHERE-2':
        event = read_electonic_event_file(name_files[files])

    #Signal with the electronic SPHERE-3
    if type_analyse=='electronic' and name_telesc=='SPHERE-3':
        event = read_electonic_event_file_sph3(name_files[files])

    if type_analyse=='experiment' and name_telesc=='SPHERE-2':
        event = read_experiment_sph2(name_files[files])
        H=float(list(pd.read_csv(name_files[files]).iloc[4])[0])-456 #430 = Baikal altitude
        b=0.9051146902370498*(H/500) 
        x_snow=-b*x_mos #x on snow
        y_snow=-b*y_mos #y on snow
        t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns

 
        
    ## START OF THE ALGORITHM

    #Select the event with full electronic response

    if event.shape[1]!=total_number_of_pmt: #Delete if there are no 109 (SPH2) or no 379 (SPH3) PMT
        continue
    if event.shape[0]!=len(cells): #Delete if there are no 1020 cells
        continue
        
   

    #_______________________

    ### SEPARATING THE SIGNAL FROM THE BACKGROUND

    #Selection criteria

    max_impulse_total=event.sum(axis=1)[1:-1].max() #Signal maximum of summed channels

    #Piedestal

    noise_total_mean=event.sum(axis=1)[1:pied][event.sum(axis=1)[1:pied] > 0].mean() #mean of piedestal
    noise_total_std=event.sum(axis=1)[1:pied][event.sum(axis=1)[1:pied] > 0].std() #std of piedestal
    noise_mean_std=noise_total_mean+3*noise_total_std #mean noise in piedestal + 3* std

    ratio_signal_noise=max_impulse_total/noise_mean_std #ratio of max signal / noise
    if type_analyse!='only_sig' and ratio_signal_noise<2:
        continue
    #_______________________

    #New time diapason of the signal

    found_diapason=diapason_impulse(event) #Using of a FUNCTION
    event_diapason,cells_diapason = found_diapason[0],found_diapason[1] #Extract the results from the function
    event_diapason.index=np.arange(len(event_diapason)) #reindex new finded diapason

    if cells_diapason[0]+30*12.5<pied*12.5:
        continue

    #_______________________

    #Lenght of the signal

    if type_analyse=='only_sig': #just the signal without background
        event_cells_sup_0=pd.DataFrame({'event':list(event.max(axis=1)),'cells':cells})
        event_cells_sup_0=event_cells_sup_0[event_cells_sup_0['event']>0]
        
        #Delete 90%
        lim_inf,lim_sup = np.percentile(event_cells_sup_0['cells'], 5), np.percentile(event_cells_sup_0['cells'], 95) #percent of the impulse saved
        idx_95=cells[event_cells_sup_0['cells'][event_cells_sup_0['cells']>lim_sup].index[0]-1]
        idx_5=cells[event_cells_sup_0['cells'][event_cells_sup_0['cells']<lim_inf].index[-1]+1]
        lenght_impulse=idx_95-idx_5 #lenght
        
    elif type_analyse!='only_sig': #signal with background
        lenght_impulse=cells_diapason[-1]-cells_diapason[0]-60*12.5  #lenght


    if lenght_impulse<100 or lenght_impulse>2000 :
        continue
    
    if max_impulse_total>10000 or lenght_impulse>2000:
        continue
    
    
    
    #_______________________


    #Amplification by sliding window method for each PMT

    #For the whole signal

    all_window=event_diapason.rolling(window=4).sum().dropna()#in the new finded diapason, N photons
    all_window.index=np.arange(0,len(all_window)) #reindex

    cells_window=cells_diapason[:-3] #for the time

    #Just for the background

    noise_window=event[1:pied].rolling(window=4).sum().dropna() #extract the piedestal
    mean_noise=noise_window[noise_window>0].mean() #mean in each PMT
    std_noise=noise_window[noise_window>0].std()  #std in each PMT

    if type_analyse=='only_sig':
        mean_std=[0]*total_number_of_pmt #threshold <bg>+3 std
    else:
        mean_std=mean_noise+3*std_noise #threshold <bg>+3 std

    #_______________________

    #Determine the axis of the shower

    #Real 

    if type_analyse!='electronic' and type_analyse!='experiment': #no information about the axis in electronic files
        x0_real=float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[2]) #real x of the axis
        y0_real=float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[3]) #real y of the axis
    if type_analyse=='electronic':
        path_elec_d0='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_30PeV_pro_500m/'+name_files[files].split('/')[-1]
        x0_real=float(pd.read_csv(path_elec_d0).columns.tolist()[0].split()[2])
        y0_real=float(pd.read_csv(path_elec_d0).columns.tolist()[0].split()[3]) #real y of the axis

    #AXIS
    
    x0_max=x_mos[event_diapason.sum().idxmax()] #on mosaic, mm
    y0_max=y_mos[event_diapason.sum().idxmax()] #on mosaic, mm
    d0_from_center=np.sqrt(x0_max**2+y0_max**2) #m
    
    if d0_from_center>d0_center:
        continue

    results_axis=axis(x0_max,y0_max,event_diapason) #Using of a FUNCTION
    x0,y0=results_axis[0]*-b,results_axis[1]*-b #return x and y of the axis in meters on snow

    d0=np.sqrt(x0**2+y0**2)

    #_______________________

    ### FIND THE SIGNAL IN EACH PMT

    #Save the PMT with the max > threshold

    df_front_N=pd.DataFrame({'PMT':num_pmt,'x':x_snow,'y':y_snow,'N_max':list(all_window.max()),'t_max':cells_window[all_window.idxmax()],'index':all_window.idxmax(),'noise_thrs':list(mean_std)}).sort_values(by='N_max',ascending=False)
    df_front_N=df_front_N[df_front_N['N_max']>df_front_N['noise_thrs']] #save the PMT up to max noise

    #Translate max of t on snow

    df_front_N['t_max']=df_front_N['t_max']-t_path[np.array(df_front_N['PMT']).astype(int)]

    #_______________________

    #Start the DFS method

    #Separate in 2 tables

    df_front_N_1=df_front_N.iloc[0].to_frame().T #Only the max from each max PMT - axis of the shower
    df_front_N_2=df_front_N.iloc[1:] #Others
     
    #Find 6 neighbors around the max

    PMT_around=circle_of_pmt_around_pmt.loc[df_front_N_1.index[0]].dropna()[1:]  #PMT around the max, delete the central PMT

    #Move them from 1 to the 2nd table

    df_front_N_1=pd.concat([df_front_N_1,df_front_N_2.loc[PMT_around[PMT_around.isin(df_front_N_2['PMT'])]]]) 
    df_front_N_2=df_front_N_2.drop(PMT_around[PMT_around.isin(df_front_N_2['PMT'])]) #Delete 6 neighbors around the max PMT

    # Return the table with saved PMT

    good_PMT=saved_PMT(df_front_N_1,df_front_N_2) #Using of a FUNCTION
    
    # good_PMT=df_front_N[df_front_N['N_max']>3*df_front_N['noise_thrs']] #save the PMT up to max noise

    # good_PMT['t_max_on_min']=np.array(good_PMT['t_max'])-min(np.array(good_PMT['t_max'])) #Readjust ti on t min
    # good_PMT=good_PMT.sort_values(by='N_max')

    good_PMT['t_max_on_min']=np.array(good_PMT['t_max'])-min(np.array(good_PMT['t_max'])) #Readjust ti on t min
    good_PMT=good_PMT.sort_values(by='N_max')

    

    # Delete PMT if there is less than 3 PMT around
    # def delete_less_3_nbr(good_PMT):
    #     bad_PMT=[]
    #     for i in range(0,len(good_PMT)):
    #         if len(good_PMT['PMT'][good_PMT['PMT'].isin(circle_of_pmt_around_pmt.loc[good_PMT.iloc[i][0]].dropna())])<3:
    #             bad_PMT.append(int(good_PMT.iloc[i][0]))
        
    #     good_PMT_new=good_PMT.drop(bad_PMT)

    #     return good_PMT_new

    good_PMT=delete_less_3_nbr(good_PMT)
    
    # good_PMT=good_PMT.sort_values(by='PMT')[0:int(len(good_PMT)/2)]

    if len(good_PMT)<10:
        continue

    # Front on the snow

    #SAVE FOR ELENA ALEKSEEVNA
    
    x_front,y_front,t_front,N_front=np.array(good_PMT['x']),np.array(good_PMT['y']),np.array(good_PMT['t_max_on_min']),np.array(good_PMT['N_max'])

    df_front=good_PMT[['PMT','x','y','t_max_on_min','N_max']].sort_values(by='PMT').rename(columns={"t_max_on_min": "t", "N_max": "N"})
    df_front['PMT']=df_front['PMT'].astype(int)
    
    filesname=name_files[files].split('/')
    filesname.insert(-1,'front')
    dirname="/".join(filesname[:-1])
    
    if type_analyse!='experiment':
        filesname="/".join(filesname)+'_front'
    else:
        filesname = "/".join(filesname).replace(".txt", "_front.txt")
        
    # if not os.path.isdir(dirname):
    #     os.mkdir(dirname)
    # df_front.to_csv(filesname, index=False)
    

    
    #END SAVE
    
    both_angles=angles(x_front,y_front,t_front,N_front)
    theta=both_angles[0]
    phi=both_angles[1]
    
    if both_angles[4]<0:
        continue
    
    # df_results=pd.DataFrame({'x_axis':x0,'y_axis':y0,'x_real':,'y_real':,'theta':theta,'phi':phi,'phi_real':,'theta_real':})
    theta_real=real_theta[int(name_files[files][-8:-5])-1]
    phi_real=real_phi[int(name_files[files][-8:-5])-1]
    if type_analyse=='experiment':
        df_results=[x0,y0,theta,phi]
    else:
        theta_real=real_theta[int(name_files[files][-8:-5])-1]
        phi_real=real_phi[int(name_files[files][-8:-5])-1]
        
        df_results=[x0,x0_real,y0,y0_real,theta,theta_real,phi,phi_real]

    
    theta_max.append(theta*180/np.pi)
    phi_max.append(phi*180/np.pi)

    theta_fit=theta
    phi_fit=phi
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
 
    delta_by_max.append(delta)

    # #____________________________________
 
    # #Return time of arrival by PART OF INTEGRAL
    
    PMT_max,index_max=np.array(good_PMT['PMT'].astype(int)),np.array(good_PMT['index'].astype(int))
    
    time_part_of_int=t_part_of_int(PMT_max,index_max)  
    t_snow_int,N_front_reconstr,x_front_int,y_front_int, PMT_int=time_part_of_int[0],time_part_of_int[1] ,time_part_of_int[2] ,time_part_of_int[3], time_part_of_int[4]
        
    #Return ANGLES by part of integral

    # t_part_int=( np.array(good_PMT['t_max_on_min'])+np.array(t_snow_int)-min(t_snow_int) ) /2
    t_part_int=np.array(t_snow_int)-min(t_snow_int)

    both_angles=angles(np.array(x_front_int),np.array(y_front_int),t_part_int,np.array(N_front_reconstr))
    
    theta=both_angles[0]
    phi=both_angles[1]
    
    df_results_all.append(df_results)
        
    # ======================================================================================
    # ======================================================================================
 
 
    # # #========================================== Delta
 
    # print('theta_fit=',theta*180/np.pi)
    # print('phi_fit=',phi*180/np.pi)
 
    # print('theta_real=',theta_real*180/np.pi)
    # print('phi_real=',phi_real*180/np.pi)
 
    theta_fit=theta
    phi_fit=phi
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
 


    files_saved.append(name_files[files])
    
    delta_by_part_of_int.append(delta)
    
    total_nbr_pmt.append(len(PMT_max))
    
    lenght_impulse_total.append(lenght_impulse)
    
    d0_all.append(d0)
 
    big_max_summed_sig.append(max_impulse_total)
    
    event_diapason_tot.append(event_diapason.sum().sum())
    
    # x0_real.append(np.nan)
    # y0_real.append(np.nan)
    
    # vall=float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[2])
    # x0_real.append(float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[2]))
    # y0_real.append(np.nan)
    
    x0_max_tot.append(x0_max)
    y0_max_tot.append(y0_max)

    x0_grav.append(x0)
    y0_grav.append(y0)

#total_d0=np.array([x0_real,y0_real,x0_max_tot,y0_max_tot,x0_grav,y0_grav]).T
# np.savetxt('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_files_elec/test',files_saved,fmt='%s')
if type_analyse=='experiment':
    axis_angles_df=pd.DataFrame(df_results_all,columns=['x_axis','y_axis','theta','phi'])
else:
    axis_angles_df=pd.DataFrame(df_results_all,columns=['x_axis','x_real','y_axis','y_real','theta','theta_real','phi','phi_real'])
###==================================================================
###==================================================================

print(np.mean(delta_by_max))
print(np.mean(delta_by_part_of_int))


# #========================== Tau

# x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
# y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
# z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
# R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
# tau=a0+a1_front*R+a2_front*R**2 +z_casc/c_ns


# ###==================================================================
# ###==================================================================

# #========= GRAPHICS

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

# ==================================

#DELTA

# plt.hist(delta_by_max)
# plt.legend(title='{} m, {} PeV, P\nby max \n$\delta$~{:.2f}º'.format(H,En,np.mean(delta_by_max)))
# plt.xlabel('$\delta$, º')
# plt.show()

plt.scatter(total_nbr_pmt,delta_by_max)
plt.ylabel('$\delta$, º')
plt.xlabel('Number of PMT')
plt.show()

plt.scatter(total_nbr_pmt,delta_by_part_of_int)
plt.ylabel('$\delta$, º')
plt.xlabel('Number of PMT')
plt.show()

plt.scatter(lenght_impulse_total,delta_by_max)
plt.ylabel('$\delta$, º')
plt.xlabel('Длина импульса')
plt.show()

plt.scatter(d0_all,delta_by_max)
plt.ylabel('$\delta$, º')
plt.xlabel('d0, m')
plt.show()


plt.scatter(big_max_summed_sig,delta_by_max)
plt.ylabel('$\delta$, º')
plt.xlabel('big max, m')
plt.show()






