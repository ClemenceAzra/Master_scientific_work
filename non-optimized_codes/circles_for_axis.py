%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit


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
En=10

###=========================================================================

##=== OPEN FILES

#===============
#== For all types of analyse

dirname = '/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/'

#PMT on mosaic: x (mm),y (mm), num
if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv(dirname + 'Initial_data/mosaic_pmt_coords_df.csv')

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv(dirname + 'Initial_data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel

#===============
#== For 1 type of analyse

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    path='/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/SPHERE-3/only_sig/500m_10PeV_P/moshits_Q2_atm01_0014_10PeV_15_000_c002'.format(En,En) #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_P_{}m/mosaic_hits_m01_pro_{}PeV_10-20_001_c001'.format(En,H,En) #Electronic signal #044 #072
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_30PeV_P_500m/mosaic_hits_m01_pro_30PeV_10-20_001_c049' #Electronic signal #044 #072
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res/mosaic_hits_m01_Fe_10PeV_10-20_038_c005' #Electronic signal #044 #072

if name_telesc=='SPHERE-3' and type_analyse=='only_sig':
    path=dirname + f'SPHERE-3/only_sig/500m_{En}PeV_P/moshits_Q2_atm01_0014_{En}PeV_15_000_c007' #01 for diploma
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_results/bad_res_sph3/moshits_Q2_atm01_0014_10PeV_15_011_c003'

###=========================================================================

##=== CONSTANT VALUES

cells=np.arange(0,1020*12.5,12.5)

if name_telesc=='SPHERE-2':
    if H==500:
        a=0.0005444285453203233 #a, b for different H --> find in another work
        b=0.9051146902370498
        nbr_pmt=109

    elif H==900:
        a=0.0008982897531551649
        b=1.6455812746636922
        
       
if name_telesc=='SPHERE-3':
    b=(1.1314175328971585)/(1000/H) #If H=500 - a --> a/2
    nbr_pmt=379

num_pmt_list=np.arange(0,nbr_pmt)



c_ns=0.3 #speed of light m/ns

part_of_int=0.5

#=====FNC

#a0 fnc
# a_a0=-0.630484 #more in the article "Алгоритм восстановления направления прихода ШАЛ для телескопа СФЕРА"
# b_a0=53.541113661

# a2=0.0005

#Paramatrization
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

def theta_initt(theta_init):
    return theta_init*len(theta_init)

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

real_theta_deg=0.1745*180/np.pi
real_phi_deg=3.0564*180/np.pi

theta_real=0.3277 #[0]
phi_real=3.0564 #[0]

# theta_real=0.3332 #[14]
# phi_real=5.7728 #[14]

# theta_real=0.2175
# phi_real=2.982

# theta_real=0.1745
# phi_real=3.0564

##=== PRELIMINARY CALCULATIONS

#==Translation x,y from mos to snow

if name_telesc=='SPHERE-2':
    x_mos=np.array(x_y_mos['x']) #x on mosaic
    y_mos=np.array(x_y_mos['y']) #y on mosaic
    x_snow=-b*x_mos #x on snow
    y_snow=-b*y_mos #y on snow

if name_telesc=='SPHERE-3':
    x_mos=np.array(pixel_data.groupby('segment').mean()['x']) #x on mosaic
    y_mos=np.array(pixel_data.groupby('segment').mean()['y']) #y on mosaic
    x_snow=-b*x_mos #x on snow
    y_snow=-b*y_mos #y on snow

#==Translation t from mos to snow (path from mos to snow)
  
t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns

#==Fonction of the approximation of the front

def tau_fnc(theta,phi,a0,a1,a2):
    x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
    y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
    z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
    R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
    return a0+a1*R+a2*R**2+z_casc/c_ns


###=========================================================================

##=== DEFINITION OF THE TYPE OF ANALYSE

#Only signal

if name_telesc=='SPHERE-2':
    num_t_event=4
if name_telesc=='SPHERE-3': 
    num_t_event=5


if type_analyse!='electronic':

    def read_mosaic_hits_file(path, step = 12.5):
        
        a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+') 
        a[num_t_event]=a[num_t_event]-min(a[num_t_event])  
        pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[num_t_event])}) #creation of DataFrame to facilite calculations
        q = np.zeros((nbr_pmt, 1020))  # empty array109 PMT * 1020 cells
        cells_here=np.arange(0,1021*step,step)
           
        shift = 500 #first bin of the signal
         
        for i in num_pmt_list:
            impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
            n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
            n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
            q[i]=n_phot_shift

#Signal with background

            if type_analyse=='sig_&_bg':
                q[i]=n_phot_shift+np.random.choice(data_photons_bg['{}'.format(i)].dropna(),size=1020)
                
        return pd.DataFrame(q.T)
    event = read_mosaic_hits_file(path)

#Electronic

if type_analyse=='electronic':

    def read_electonic_event_file(path):
        event = pd.read_csv(path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip')
        for pmt in range(event.shape[1]):
            event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
            event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
        return event
    event = read_electonic_event_file(path)
    
H=500
if name_telesc=='SPHERE-2':
    
    #Coeff for translate x mos --> x snow
    b=0.9051146902370498*(H/500) 
    
    #Values for the axis
    # d0_center=200*(H/500)  #last circle on mosaic
    # d0_before_c=155*(H/500) #before last circle on mosaic
    # d_six=50*(H/500)  #distance PMT and 1 circle around
    # d_double_six=100*(H/500) #distance PMT and 2 circles around
    d0_center=192 #last circle on mosaic
    d0_before_c=154 #before last circle on mosaic
    d_six=45  #distance PMT and 1 circle around
    d_double_six=95 #distance PMT and 2 circles around
    
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
###========================================================================


#Search the impulse in the отклик

pulse = event.sum(axis=1) # total impulse in 1020 cells
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
cells_impulse_filtrate_each_pmt=[]
impulse_filtrate_each_pmt=[]

for i in range(0,event_diapason.shape[1]):
    event_pmt=event_diapason[i]
    
    #Condition 1 for PMT: delete the PMT with max signal < max BG
    
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

    #Condition 2 for bins: delete bins < max BG
    
    impulse_filtrate=impulse[impulse>max(event[i][0:400])] 
    cells_impulse_filtrate=cells_impulse[impulse>max(event[i][0:400])] 

        
   #Condition 3 for bins: delete the PMT with 0 bins
   
    if len(cells_impulse_filtrate)<1: #if number of bins in PMT less than 4 - low signal
        continue
    
    #Aoment of the arrival of shower in the PMT
    t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

    #SAVE DATA
    
    t_front_reconstr.append(t_arrival-t_path[i]) #t - obtain front
    x_front_reconstr.append(x_snow[i]) #x - obtain front
    y_front_reconstr.append(y_snow[i]) #y - obtain front
    
    saved_pmt.append(num_pmt_list[i]) # Nº of the saved PMT with signal
    number_photons_i.append(impulse_filtrate.sum()) #number of signal photons in each PMT

    cells_impulse_filtrate_each_pmt.append(cells_impulse_filtrate) #to see image
    impulse_filtrate_each_pmt.append(impulse_filtrate) #to see image

###============================================

#Axis of the shower

# number_photons=np.array(number_photons)


x0_front=x_front_reconstr[number_photons_i.index(max(number_photons_i))] 
y0_front=y_front_reconstr[number_photons_i.index(max(number_photons_i))]
t0_front=t_front_reconstr[number_photons_i.index(max(number_photons_i))]

d0_from_center=np.sqrt(x0_front**2+y0_front**2)

d_x_to_x_max=np.sqrt( (np.array(x_front_reconstr)-x0_front)**2  + (np.array(y_front_reconstr)-y0_front)**2  )

# if H==500:
#     if d0_from_center>140 and d0_from_center<180:
#         d_max=50
#     if d0_from_center<140:
#        d_max=100

       
# if H==900:
#     if d0_from_center>2 and d0_from_center<330:
#         d_max=100
#     if d0_from_center<260:
#        d_max=150

# if d0_from_center>144*(H/500) and d0_from_center<d0_center*(H/500):
#     d_max=30*(H/500)
# if d0_from_center<120*(H/500):
#    d_max=100*(H/500)   

# x_circle_around_max=np.array(x_front_reconstr)[np.where(d_x_to_x_max<d_max)[0]]
# y_circle_around_max=np.array(y_front_reconstr)[np.where(d_x_to_x_max<d_max)[0]]
# number_photons_around_max=np.array(number_photons_i)[np.where(d_x_to_x_max<d_max)[0]]

# x0_grav_center=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
# y0_grav_center=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)

# x0_real=float(pd.read_csv(path).columns.tolist()[0].split()[2])
# y0_real=float(pd.read_csv(path).columns.tolist()[0].split()[3])

# all_window=event_diapason.rolling(window=4).sum().dropna()
# all_window.index=np.arange(0,len(all_window))

# x0_center=x_mos[all_window.max().idxmax()]
# y0_center=y_mos[all_window.max().idxmax()]

# print(x0_center,y0_center)

###==================================================================

#========= GRAPHICS

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

#Mosaic with last circle, before last circle, and interior

# plt.scatter(x_y_mos['x'],x_y_mos['y'])
# plt.scatter(-x0_center,-y0_center)
# plt.axis('square')
# plt.show()

#==========

#3 ZONES OF MOSAIC

# exterior=d0_center*(H/500)
# before_last=d0_before_c*(H/500)

d_x_to_center=np.sqrt( np.array(x_mos)**2  + np.array(y_mos)**2  )
 
x_last_circle=np.array(x_mos)[np.where(d_x_to_center>d0_center)[0]]
y_last_circle=np.array(y_mos)[np.where(d_x_to_center>d0_center)[0]]

x_before_last_circle=np.array(x_mos)[np.where(d_x_to_center>d0_before_c)[0]]
y_before_last_circle=np.array(y_mos)[np.where(d_x_to_center>d0_before_c)[0]]

plt.scatter(x_mos,y_mos,c='greenyellow',label='interior')
plt.scatter(x_before_last_circle,y_before_last_circle,label='before \nlast circle',color='blue')
plt.scatter(x_last_circle,y_last_circle,label='last circle',color='orange')
plt.axis('square')
plt.xlim(-350,350)
plt.ylim(-350,350)
plt.xlabel('x на мозаике, мм')
plt.ylabel('y на мозаике, мм')
plt.legend(fontsize=8)
plt.show()

#Circle around

#Modelisation of 3 max

if name_telesc=='SPHERE-2':
    #Snow
    x0_1, y0_1=-20.311293184751595, -35.180191474677365 #interior
    x0_2, y0_2=-79.79018458238589, 138.2006496359943 #before last
    x0_3, y0_3=39.903004803815, 204.21223889454336 #last
    
    # mosaic
    x0_1, y0_1=-22.186325, -115.76552 #interior
    x0_2, y0_2=-88.154778,152.68855 #before last
    x0_3, y0_3= 44.086131, 225.62029 #last
    
    
if name_telesc=='SPHERE-3':
    #Snow
    x0_1, y0_1=-20.311293184751595, -35.180191474677365 #interior
    x0_2, y0_2=-79.79018458238589, 138.2006496359943 #before last
    x0_3, y0_3=39.903004803815, 204.21223889454336 #last
    
    # mosaic
    x0_1, y0_1=79.58326271428572, -31.808191285714287 #interior
    x0_2, y0_2=-282.80600999999996,-10.363606997142856 #before last
    x0_3, y0_3= 23.790305714285715, -311.2207528571429 #last


d_x_to_x_max_1=np.sqrt( (np.array(x_mos)-x0_1)**2  + (np.array(y_mos)-y0_1)**2  )
d_x_to_x_max_2=np.sqrt( (np.array(x_mos)-x0_2)**2  + (np.array(y_mos)-y0_2)**2  )

x_circle_around_max_1=np.array(x_mos)[np.where(d_x_to_x_max_1<d_double_six)[0]] #interior
y_circle_around_max_1=np.array(y_mos)[np.where(d_x_to_x_max_1<d_double_six)[0]] #interior

x_circle_around_max_2=np.array(x_mos)[np.where(d_x_to_x_max_2<d_six)[0]] #before last
y_circle_around_max_2=np.array(y_mos)[np.where(d_x_to_x_max_2<d_six)[0]] #before last

# 2 types of circle around 2 different max

plt.scatter(x_mos,y_mos)
plt.scatter(x_circle_around_max_1,y_circle_around_max_1,c='red')
plt.scatter(x_circle_around_max_2,y_circle_around_max_2,c='red')
plt.scatter(x0_1,y0_1,c='lime')
plt.scatter(x0_2,y0_2,c='darkblue')
plt.scatter(x0_3,y0_3,c='orange')
plt.axis('square')
plt.xlim(-350,350)
plt.ylim(-350,350)
plt.xlabel('x на снегу, мм')
plt.ylabel('y на снегу, мм')
plt.legend(fontsize=9)
plt.show()


# d0_center=d0_center*(H/500)
# before_last=d0_before_c*(H/500)

d_x_to_center=np.sqrt( np.array(x_mos)**2  + np.array(y_mos)**2  )
 
x_last_circle=np.array(x_mos)[np.where(d_x_to_center>d0_center)[0]]
y_last_circle=np.array(y_mos)[np.where(d_x_to_center>d0_center)[0]]

x_before_last_circle=np.array(x_mos)[np.where(d_x_to_center>d0_before_c)[0]]
y_before_last_circle=np.array(y_mos)[np.where(d_x_to_center>d0_before_c)[0]]

marker_ok="o"
m = MarkerStyle(marker_ok)
m._transform.rotate_deg(75)

# plt.scatter(x_mos,y_mos,c='greenyellow',label='Внутренняя \nчасть')
# plt.scatter(x_before_last_circle,y_before_last_circle,label='Предпоследнее \nкольцо',color='blue')
# plt.scatter(x_last_circle,y_last_circle,label='Последнее \nкольцо',color='orange')

plt.scatter(x_mos,y_mos,linewidth=2.2,marker=m,s=size_dots,edgecolor = 'dimgray',color='white')

# plt.scatter(x_before_last_circle,y_before_last_circle,linewidth=2.2,marker=m,s=size_dots,edgecolor = 'dimgray',color='white')
# plt.scatter(x_last_circle,y_last_circle,linewidth=2.2,marker=m,s=size_dots,edgecolor = 'dimgray',color='white')

plt.scatter(x_mos,y_mos,c='greenyellow',label='Interior \npart') #greenyellow
plt.scatter(x_before_last_circle,y_before_last_circle,label='Penultimate \nring',color='orange')
plt.scatter(x_last_circle,y_last_circle,label='Last \nring',color='darkred')


plt.scatter(x_circle_around_max_1,y_circle_around_max_1,c='blue')
plt.scatter(x_circle_around_max_2,y_circle_around_max_2,c='blue')
plt.scatter(x0_1,y0_1,c='lime')
plt.scatter(x0_2,y0_2,c='orange')
plt.scatter(x0_3,y0_3,c='darkred')
plt.axis('square')
plt.xlim(-400,400)
plt.ylim(-400,400)
# plt.xlabel('x, мм',fontsize=14)
# plt.ylabel('y, мм',fontsize=14)
plt.xlabel('x, mm',fontsize=20)
plt.ylabel('y, mm',fontsize=20)
plt.legend(fontsize=14, loc='center left', bbox_to_anchor=(1, 0.5))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20, ticks = np.arange(-400, 400+200, 200))
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/POS/mosaic_for_axis_sph2.eps',bbox_inches='tight')
plt.show()

#500 m 30 PeV for SPHERE-2, 1000 m 10 PeV for SPHERE-3


# plt.scatter(x_mos,y_mos,c='greenyellow',label='Внутренняя \nчасть')
# plt.scatter(x0_center,y0_center)
# plt.axis('square')
# plt.show()






