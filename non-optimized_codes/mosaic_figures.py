import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


###=========================================================================
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
###=========================================================================

dirname = '/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/Initial_data/'

if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv(dirname + 'mosaic_pmt_coords_df.csv')
    PMT=list(x_y_mos['pmt'])

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv(dirname + 'SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel
    x_y_mos=pd.DataFrame({'x':pixel_data.groupby('segment').mean()['x'],'y':pixel_data.groupby('segment').mean()['y']})
    PMT=list(pixel_data.groupby('segment').mean().index)
###=========================================================================
###=========================================================================

if name_telesc=='SPHERE-2':
    if H==500:
        a=0.0005444285453203233 #a, b for different H --> find in another work
        b=0.9051146902370498
    elif H==900:
        a=0.0008982897531551649
        b=1.6455812746636922
       
if name_telesc=='SPHERE-3':
    b=(1.1314175328971585)/(1000/H) #If H=500 - a --> a/2

###=========================================================================
###=========================================================================

if name_telesc=='SPHERE-2':
    x_snow=-b*x_y_mos['x']
    y_snow=-b*x_y_mos['y']
    size_dots=80
    lim_size=400
    size_number=3
    marker_ok="o"

if name_telesc=='SPHERE-3':
    x_snow=-b*np.array(pixel_data.groupby('segment').mean()['x'])
    y_snow=-b*np.array(pixel_data.groupby('segment').mean()['y'])
    size_dots=95
    size_number=2.5
    lim_size=400
    marker_ok="h"

distance=np.sqrt( np.array(x_y_mos['x'])**2 + np.array(x_y_mos['y'])**2)
    
###=========================================================================
###=========================================================================

# #========= GRAPHICS
from matplotlib.markers import MarkerStyle

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

#ON MOSAIC
m = MarkerStyle(marker_ok)
m._transform.rotate_deg(75)

plt.scatter(list(x_y_mos['x']),list(x_y_mos['y']),linewidth=2.2,marker=m,s=size_dots,edgecolor = 'dimgray',color='white')

plt.axis('scaled')
plt.xlabel('mm',fontsize=18)
plt.ylabel('mm',fontsize=18)
for x,y,z in zip(list(x_y_mos['x']),list(x_y_mos['y']),PMT):
    plt.annotate(int(z), # this is the text
                  (x,y), # these are the coordinates to position the label
                  textcoords="offset points", # how to position the text
                  xytext=(0,-1.5), # distance from text to points (x,y)
                  ha='center',color='darkred',size=size_number)
plt.ylim(-lim_size,lim_size)
plt.xlim(-lim_size,lim_size)
plt.xticks(fontsize=18, ticks = np.arange(-400, 400+200, 200))
plt.yticks(fontsize=18, ticks = np.arange(-400, 400+200, 200))
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Постер ECCR/Figures/mosaic_sph2.pdf',bbox_inches='tight')
plt.show()


#ON SNOW

# plt.scatter(x_snow,y_snow,marker="h",s=450,edgecolor = 'black',color='white')
# plt.axis('equal')
# plt.xlabel('x on snow, m')
# plt.ylabel('y on snow, m')
# for x,y,z in zip(x_snow,y_snow,PMT):
#     plt.annotate(int(z), # this is the text
#                  (x,y), # these are the coordinates to position the label
#                  textcoords="offset points", # how to position the text
#                  xytext=(0,-3), # distance from text to points (x,y)
#                  ha='center',color='darkred',size=8)
# plt.legend(title='H={} m'.format(H))
# plt.show()

