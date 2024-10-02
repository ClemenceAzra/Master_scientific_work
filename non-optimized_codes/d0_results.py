import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

###=============================================================
###=============================================================

## !!! CHANGE ONLY HERE

# nuclei='Fe'
# nuclei='N'
nuclei='P'

# Choose the altitude
H=900 #m, 900 

# Choose the energy
En=10 #PeV, 10

integral=0.5

# Choose the analyse
analyse='only_sig' #elec, bg
# analyse='electronic' #elec, bg
telescope='sph2'
angle_th='axis_pmt'

dirname='/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/'

## !!! END CHANGE ONLY HERE

###=============================================================
###=============================================================

#==================================================================
if telescope!='sph3':
    all_val_x0_y0=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/sph2/d0/{}/{}_d0_{}m_{}PeV_{}_sph2_{}.csv'.format(analyse,analyse,H,En,nuclei,angle_th),header=None,sep='\s+')
    vision=180
    in_legend='СФЕРА-2'
    if H==500:
        vision=180
        max_x_axis=180
    elif H==900:
        vision=310
        max_x_axis=330

if telescope=='sph3':
    all_val_x0_y0=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/{}/d0/{}/{}_d0_{}m_{}PeV_{}_sph3_{}.csv'.format(telescope,analyse,analyse,H,En,nuclei,angle_th),header=None,sep='\s+')
    in_legend='СФЕРА-3'
    if H==500:
        vision=165
        max_x_axis=180
    elif H==1000:
        vision=330
        max_x_axis=330



if analyse=='electronic':
    
    files=list(pd.read_csv(dirname+'RESULTS/ALL_RESULTS/files/{}/{}_files_{}m_{}PeV_{}.csv'.format(analyse,analyse,H,En,nuclei),header=None)[0])
    order_elec=[ int(files[i][-8:-5] + files[i][-3:]) for i in range(0,len(files))]
    all_val_x0_y0['order_elec']=order_elec
    all_val_x0_y0=all_val_x0_y0.sort_values(by=['order_elec'])

    real_d0=pd.read_csv(dirname+'/RESULTS/ALL_RESULTS/d0/d0_real_elec/d0_real_{}m_{}PeV_{}.csv'.format(H,En,nuclei),header=None,sep='\s+')
    real_d0=real_d0.sort_values(by=[2])

    all_val_x0_y0[0]=real_d0[0]
    all_val_x0_y0[1]=real_d0[1]
    
    in_legend='СФЕРА-2'
    
    
if analyse=='only_sig' and  telescope!='sph2':
    
    all_val_x0_y0=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/{}/d0/{}/{}_d0_{}m_{}PeV_{}_{}_{}.csv'.format(telescope,analyse,analyse,H,En,nuclei,telescope,angle_th),header=None,sep='\s+')
    
    in_legend='SPHERE-3'
    
    
    
# elif analyse=='only_sig' and telescope!='sph3':
#     in_legend='Чистый сигнал'

    
#Extract values
  
x0_real=np.array(all_val_x0_y0[0])
y0_real=np.array(all_val_x0_y0[1])

x0_max=np.array(all_val_x0_y0[2])
y0_max=np.array(all_val_x0_y0[3])

x0_grav_center=np.array(all_val_x0_y0[4]) 
y0_grav_center=np.array(all_val_x0_y0[5]) 

d0_real=np.sqrt(x0_real**2+y0_real**2)
d0_grav_center=np.sqrt(x0_grav_center**2+y0_grav_center**2)

#ERROR


def err_d0(results):
    delta_axis_total=[]
    d0_real=np.sqrt(x0_real**2+y0_real**2)
    d0_max=np.sqrt(x0_max**2+y0_max**2)
    d0_grav_center=np.sqrt(x0_grav_center**2+y0_grav_center**2)

    delta_d0_max=abs(d0_real-d0_max)
    delta_d0_grav_center=abs(d0_real-d0_grav_center)
    
    d0=pd.DataFrame({'d0_real':d0_real,'d0_max':d0_max,'d0_grav_center':d0_grav_center,'delta_max':delta_d0_max,'delta_grav':delta_d0_grav_center}).sort_values(by='d0_real')
    lenght=np.arange(0,vision+10,10)
    delta_axis_total.append(delta_d0_max)
    x_d0_grav=[]
    x_d0_max=[]
    delta_grav_mean=[]
    delta_grav_std=[]
    delta_max_mean=[]
    delta_max_std=[]
    for i in range(0,len(lenght)-1):
        delta_max=list(d0['delta_max'][ (d0['d0_real']>lenght[i]) & (d0['d0_real']<lenght[i+1]) ])
        delta_grav=list(d0['delta_grav'][ (d0['d0_real']>lenght[i]) & (d0['d0_real']<lenght[i+1]) ])
        if len(delta_grav)>2:
            x_d0_grav.append(lenght[i])
            delta_grav_mean.append(np.mean(delta_grav))
            delta_grav_std.append(np.std(delta_grav))
        if len(delta_max)>2:
            x_d0_max.append(lenght[i])    
            delta_max_mean.append(np.mean(delta_max))
            delta_max_std.append(np.std(delta_max))
            

    return delta_grav_mean, delta_max_mean, delta_grav_std,delta_max_std,x_d0_grav,x_d0_max, delta_axis_total

#============

results_d0=err_d0(all_val_x0_y0)
delta_grav_mean,delta_max_mean,delta_grav_std,delta_max_std,x_d0_grav,x_d0_max=results_d0[0],results_d0[1],results_d0[2],results_d0[3],results_d0[4],results_d0[5]

#Normalize capsize to min =0

verify_if_below_0_max=np.array(delta_max_mean)-np.array(delta_max_std)
verify_if_below_0_grav=np.array(delta_grav_mean)-np.array(delta_grav_std)

df_verify_below_0=pd.DataFrame({'delta_max_std':delta_max_std,'verify_if_below_0_max':verify_if_below_0_max,'delta_grav_std':delta_grav_std,'verify_if_below_0_grav':verify_if_below_0_grav})

#=====================================================================

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams["font.family"] = "Times New Roman"


# plt.scatter(x_d0_max,delta_max_mean,color='red',label='По макс. ФЭУ  <{:.0f}>±{:.0f}м'.format(np.mean(delta_max_mean),np.mean(delta_max_std)), edgecolor ='black',s=15)
# plt.scatter(x_d0_grav,delta_grav_mean,color='green',label='По центру тяжести <{:.0f}>±{:.0f}м'.format(np.mean(delta_grav_mean),np.mean(delta_grav_std)), edgecolor ='black',s=15)

# plt.errorbar(x_d0_max,delta_max_mean, yerr = delta_max_std,fmt ='o',markersize=3, capsize=4,color='red')
# plt.errorbar(x_d0_grav,delta_grav_mean, yerr = delta_grav_std,fmt ='o',markersize=3, capsize=4,color='green')
# plt.ylabel(' < |∆d0| >, м',fontsize=16)
# plt.xlabel('real d0, м',fontsize=16)
# plt.legend(title='{}, {}ПэВ, {} м, {}'.format(nuclei,En,H,in_legend), loc='upper center',fontsize=12,title_fontsize=14)
# plt.xlim(-5,max_x_axis)
# plt.ylim(0,60)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.show()

#_________________________________

#d0 only by gravity center for one type of telescope
 
val=3

# plt.scatter(x_d0_grav,delta_grav_mean,color='orange',label='{} m {} PeV,  {}={:.0f}±{:.0f}m'.format(H,En,'$∆d_{tot}$',np.mean(delta_grav_mean),np.mean(delta_grav_std)), edgecolor ='black',s=15)
# plt.errorbar(x_d0_grav,delta_grav_mean, yerr = delta_grav_std,fmt ='o',markersize=1, capsize=3,color='orange')

# save_data=pd.DataFrame({'500m_x':x_d0_grav_500+[np.nan]*15,'500m_mean':delta_grav_mean_500+[np.nan]*15,'500m_std':delta_grav_std_500+[np.nan]*15,'1000m_x':x_d0_grav,'1000m_mean':delta_grav_mean,'1000m_std':delta_grav_std})
# save_data.to_csv('/Users/clemence/Documents/data_axis.csv')

# plt.scatter(np.array(x_d0_grav_500)+val,delta_grav_mean_500,color='green',label='500 m 30 PeV,  {}={:.0f}±{:.0f}m'.format('$∆d_{tot}$',np.mean(delta_grav_mean_500),np.mean(delta_grav_std_500)), edgecolor ='black',s=15)
# plt.errorbar(np.array(x_d0_grav_500)+val,delta_grav_mean_500, yerr = delta_grav_std_500,fmt ='o',markersize=1, capsize=3,color='green')


# x_d0_grav_500,delta_grav_mean_500,delta_grav_std_500=x_d0_grav,delta_grav_mean,delta_grav_std

# plt.scatter(np.array(x_d0_grav_500)+val,delta_grav_mean_500,color='green',label='  500 м, 30 ПэВ, {}={:.0f}±{:.0f}м'.format('$∆d_{tot}$',np.mean(delta_grav_mean_500),np.mean(delta_grav_std_500)), edgecolor ='black',s=15)
# plt.errorbar(np.array(x_d0_grav_500)+val,delta_grav_mean_500, yerr = delta_grav_std_500,fmt ='o',markersize=1, capsize=3,color='green')
   
plt.scatter(x_d0_grav,delta_grav_mean,color='orange',label='{} м, 10 ПэВ, {}={:.0f}±{:.0f}м'.format(H,'$∆d_{tot}$',np.mean(delta_grav_mean),np.mean(delta_grav_std)), edgecolor ='black',s=15)
plt.errorbar(x_d0_grav,delta_grav_mean, yerr = delta_grav_std,fmt ='o',markersize=1, capsize=3,color='orange')

# save_data=pd.DataFrame({'500m_x':x_d0_grav_500+[np.nan]*12,'500m_mean':delta_grav_mean_500+[np.nan]*12,'500m_std':delta_grav_std_500+[np.nan]*12,'900m_x':x_d0_grav,'900m_mean':delta_grav_mean,'900m_std':delta_grav_std})
# save_data.to_csv('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Известия_РАН/axis_sh2.csv')


plt.ylabel('∆d, м',fontsize=15)
plt.xlabel('d, м',fontsize=15)
plt.legend(title='p, СФЕРА-2', loc='upper center',fontsize=14,title_fontsize=14)
plt.ylim(0,36.5)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.rcParams['axes.axisbelow'] = True
# plt.grid(True)
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Статья_3_УЗ/Figs/d0_sphere2.eps',bbox_inches='tight')
plt.show()


# #To Elena Alekseevna

# axis_results=pd.read_csv('/Users/clemence/Documents/Магистратура_наука/data_axis.csv')

# plt.scatter(axis_results['1000m_x'],axis_results['1000m_mean'],color='orange',label='1000 м 10 ПэВ,  {}={:.0f}±{:.0f}m'.format('$∆d_{tot}$',axis_results['1000m_mean'].mean(),axis_results['1000m_std'].mean()), edgecolor ='black',s=15)
# plt.errorbar(axis_results['1000m_x'],axis_results['1000m_mean'], yerr = axis_results['1000m_std'],fmt ='o',markersize=1, capsize=3,color='orange')

# plt.scatter(axis_results['500m_x'],axis_results['500m_mean'],color='green',label='500 м 10 ПэВ,  {}={:.0f}±{:.0f}m'.format('$∆d_{tot}$',axis_results['500m_mean'].mean(),axis_results['500m_std'].mean()), edgecolor ='black',s=15)
# plt.errorbar(axis_results['500m_x'],axis_results['500m_mean'], yerr = axis_results['500m_std'],fmt ='o',markersize=1, capsize=3,color='green')

# plt.ylabel('∆d, м',fontsize=15)
# plt.xlabel('d, м',fontsize=15)
# plt.legend(title='p, СФЕРА-3', loc='upper center',fontsize=14,title_fontsize=14)
# plt.ylim(0,36.5)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# # plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Статья_3_УЗ/Figs/d0_sphere3.eps',bbox_inches='tight')
# plt.show()



