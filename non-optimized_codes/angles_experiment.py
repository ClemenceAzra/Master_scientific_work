import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statistics

year='both_years' #both_years
# year='2013' #both_years
dirname='/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/'

###=============================================================
###=============================================================

#Obtain results

results=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_total_{}_2_100.csv'.format(year),header=None,sep='\s+')
# results=results[results[4]>20]
files=list(pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_files_{}_2_100.csv'.format(year),header=None)[0])
H=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_H_{}_2_100.csv'.format(year),header=None,sep='\s+')
# H=H[results[4]>20]

#====================

#ALL 2013 and 2012

H_2013=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_H_2013_2_100.csv'.format(year),header=None,sep='\s+')
H_2012=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_H_2012_2_100.csv'.format(year),header=None,sep='\s+')

results_2012=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_total_2012_2_100.csv'.format(year),header=None,sep='\s+')
results_2013=pd.read_csv(dirname+'RESULTS/ALL_RESULTS/experiment/experiment_total_2013_2_100.csv'.format(year),header=None,sep='\s+')
d0_2012=results_2012[7]
d0_2013=results_2013[7]
length_2012=results_2012[5]
length_2013=results_2013[5]
ratio_or_max_2012=results_2012[8]
ratio_or_max_2013=results_2013[8]

#====================


theta=results[0]*180/np.pi
phi=results[1]*180/np.pi
 
nmbr_pmt=results[4]
length=results[5]
num_files=results[6]
d0=results[7]
ratio_or_max=results[8]
a0,a1,a2=results[9],results[10],results[11]

pas_theta=5
pas_phi=30
pas_length=100
pas_axis=25
pas_H=50
pas_ratio=5

angles_theta=np.arange(0,70+pas_theta,pas_theta)
angles_phi=np.arange(0,360+2*pas_phi,pas_phi)
lenght_imp=np.arange(200,1500+pas_length,pas_length)
d0_axis=np.arange(0,350+pas_axis,pas_axis)
H_high=np.arange(200,1000+pas_H,pas_H)
ratio_bg=np.arange(0,50+pas_ratio,pas_ratio)

###=============================================================
###=============================================================


def angles(angle,angles_bar):
    angles_perc=[]
    for i in range(0,len(angles_bar)-1):
        # delta_df=pd.DataFrame({'delta':delta_max})
        angles_perc.append(len(angle[ (angle>angles_bar[i]) & (angle<angles_bar[i+1])])*100/len(angle))
    return angles_perc

#Dependance front param (d0)

def analyse_param(results,x_analyse,y_analyse,step):
    colonne_d0=x_analyse
    table=results[[y_analyse,colonne_d0]] #dataframe of param a (col 1) and d0 (col 2)
    lenght=np.arange(results[colonne_d0].min(),results[colonne_d0].max()+step,step)
    x_step=[]
    y_mean=[]
    y_std=[]
    for i in range(0,len(lenght)-1):
        y=table[ (table[colonne_d0]>lenght[i]) & (table[colonne_d0]<lenght[i+1]) ][y_analyse]
        if len(y)>2:
            x_step.append(lenght[i])
            y_mean.append(y.mean())
            y_std.append(y.std())
    return  x_step, y_mean, y_std


###=============================================================
###=============================================================

results_angles_th=angles(theta,angles_theta)
results_angles_ph=angles(phi,angles_phi)

results_length=angles(length,lenght_imp)
results_length_2012=angles(length_2012,lenght_imp)
results_length_2013=angles(length_2013,lenght_imp)

results_axis=angles(d0,d0_axis)
results_axis_2012=angles(d0_2012,d0_axis)
results_axis_2013=angles(d0_2013,d0_axis)


results_ratio=angles(ratio_or_max,ratio_bg)
results_ratio_2012=angles(ratio_or_max_2012,ratio_bg)
results_ratio_2013=angles(ratio_or_max_2013,ratio_bg)

results_H=angles(H[0],H_high)
results_H_2013=angles(H_2013[0],H_high)
results_H_2012=angles(H_2012[0],H_high)

#________

#A0, A1, A2

#Dependance:(a0)
pas_param_front=50
delta_a0=analyse_param(results,7,9,pas_param_front) #total table, x, y, step
x_a0,y_a0,y_a0_std=delta_a0[0],delta_a0[1],delta_a0[2]

#Dependance:(a1)

delta_a1=analyse_param(results,7,10,pas_param_front)
x_a1,y_a1,y_a1_std=delta_a1[0],delta_a1[1],delta_a1[2]

#Dependance:(a2)

delta_a2=analyse_param(results,7,11,pas_param_front)
x_a2,y_a2,y_a2_std=delta_a2[0],delta_a2[1],delta_a2[2]

#________

#theta (length)

theta_length=analyse_param(results,5,0,100) 
x_th_length,y_th_length,y_th_length_std=theta_length[0],theta_length[1],theta_length[2]

# length (nbr pmt)

length_pmt=analyse_param(results,5,4,100)
x_length_pmt,y_length_pmt,y_length_pmt_std=length_pmt[0],length_pmt[1],length_pmt[2]

# theta (ratio)

theta_ratio=analyse_param(results,4,0,10)
x_theta_ratio,y_theta_ratio,y_theta_ratio_std=theta_ratio[0],theta_ratio[1],theta_ratio[2]


###===

H_show=len(H[  (H>300)& (H<400) ].dropna())*100/len(results)
d0_show=len(d0[  (d0>100)& (d0<200) ].dropna())*100/len(results)
ratio_show=len(d0[  (ratio_or_max>0)& (ratio_or_max<10) ].dropna())*100/len(results)

###=============================================================
###=============================================================

##==== GRAPHICS

params = {"xtick.direction": "in", "ytick.direction": "in"} #ticks in box
plt.rcParams.update(params)


# plt.hist(theta)
# plt.xlabel('$\Theta$,º')
# plt.show()


# plt.hist(phi)
# plt.xlabel('$\phi$,º')
# plt.show()

# plt.hist(lenght)
# plt.xlabel('Длина импульса, нс')
# plt.show()



#THETA AND PHI

# plt.bar(angles_theta[:-1]+pas_theta/2,results_angles_th,color='blue',width=pas_theta,alpha=0.3)
# plt.plot(angles_theta[:-1],results_angles_th,color='blue',linewidth=2,drawstyle="steps-post")

x_values=list(np.arange(0,60+5,5))
plt.hist(list(theta),color='dodgerblue',alpha=0.3,bins=x_values, edgecolor='blue')
plt.plot(x_values,list(np.histogram(theta,bins=x_values)[0])+[0],color='blue',linewidth=2,drawstyle="steps-post")

plt.xlabel('$\Theta$,º',fontsize=18)
plt.ylabel('N',fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(0,60)
plt.rcParams['axes.axisbelow'] = True
plt.grid(True)
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Постер ECCR/Figures/exp_theta.pdf',bbox_inches='tight')
plt.show()


# plt.bar(angles_phi[:-1]+pas_phi/2,results_angles_ph,color='blue',width=pas_phi,alpha=0.3)
# plt.plot(angles_phi[:-1],results_angles_ph,color='blue',linewidth=2,drawstyle="steps-post")
x_values=list(np.arange(0,360+60,60))
plt.hist(list(phi),color='dodgerblue',alpha=0.3,bins=x_values, edgecolor='blue')
plt.plot(x_values,list(np.histogram(phi,bins=x_values)[0])+[0],color='blue',linewidth=2,drawstyle="steps-post")
# plt.errorbar(np.array(x_values)+15 ,list(np.histogram(phi,bins=x_values)[0])+[0], yerr=np.sqrt(np.array(list(np.histogram(phi,bins=x_values)[0])+[0])), fmt='none',capsize=5, color='black')
plt.xlabel('$\phi$,º',fontsize=18)
plt.ylabel('N',fontsize=18)
plt.xticks(fontsize=18, ticks = np.arange(0, 360+60, 60))
plt.yticks(fontsize=18, ticks = [15, 30, 45, 60, 75])
plt.xlim(0,360)
plt.rcParams['axes.axisbelow'] = True
plt.grid(True)
plt.ylim(0, None)
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Постер ECCR/Figures/exp_phi.pdf',bbox_inches='tight')
plt.show()


# ## == LENGTH OF THE IMPULSE

# # lenght_imp_2013,results_length_2013=lenght_imp,results_length
# # plt.plot(lenght_imp_2013[:-1],results_length_2013,color='navy',linewidth=2,drawstyle="steps-post",linestyle='-',label='2013')

# plt.plot(lenght_imp[:-1],results_length_2012,color='navy',linewidth=2,drawstyle="steps-post",linestyle='-.',label='2012')
# plt.plot(lenght_imp[:-1],results_length_2013,color='navy',linewidth=2,drawstyle="steps-post",linestyle='-',label='2013')
# plt.xlabel('Длина импульса, нс',fontsize=16)
# plt.ylabel('Доля событий, %',fontsize=16)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.legend(title='Год',title_fontsize=14,fontsize=14)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_length.pdf')
# plt.show()

# # == AXIS

# d0_axis_2013,results_axis_2013=d0_axis,results_axis
# plt.plot(d0_axis_2013[:-1],results_axis_2013,color='darkorange',linewidth=2,drawstyle="steps-post",linestyle='-',label='2013')

# plt.plot(d0_axis[:-1],results_axis_2012,color='darkorange',linewidth=2,drawstyle="steps-post",linestyle='-.',label='2012')
# plt.plot(d0_axis[:-1],results_axis_2013,linewidth=2,color='darkorange',drawstyle="steps-post",linestyle='-',label='2013')

# plt.xlabel('d, м',fontsize=16)
# plt.ylabel('Доля событий, %',fontsize=16)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.legend(title='Год',title_fontsize=14,fontsize=14)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_axis.pdf')
# plt.show()



# # # # # == ALTITUDE OF THE DETECTOR, H

# # # # H_high_2013,results_H_2013=H_high,results_H
# # # # plt.plot(H_high_2013[:-1],results_H_2013,color='black',linewidth=2,drawstyle="steps-post",linestyle='-',label='2013')

# plt.plot(H_high[:-1]+pas_H/2,results_H_2012,linewidth=2,drawstyle="steps-post",color='black',label='2012',linestyle='-.')
# plt.plot(H_high[:-1]+pas_H/2,results_H_2013,linewidth=2,drawstyle="steps-post",color='black',label='2013',linestyle='-')

# plt.xlabel('H, м',fontsize=16)
# plt.ylabel('Доля событий, %',fontsize=16)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.legend(title='Год',title_fontsize=14,fontsize=14)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_H.pdf')
# plt.show()


# # # == RATIO MAX / BG

# # ratio_bg_2013,results_ratio_2013=ratio_bg,results_ratio
# # plt.plot(ratio_bg_2013[:-1],results_ratio_2013,color='green',linewidth=2,drawstyle="steps-post",linestyle='-',label='2013')

# plt.plot(ratio_bg[:-1],results_ratio_2012,color='green',linewidth=2,drawstyle="steps-post",linestyle='-.',label='2012')
# plt.plot(ratio_bg[:-1],results_ratio_2013,color='green',linewidth=2,drawstyle="steps-post",linestyle='-',label='2013')
# plt.xlabel('ratio',fontsize=16)
# plt.ylabel('Доля событий, %',fontsize=16)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.legend(title='Год',title_fontsize=14,fontsize=14)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_ratio.pdf')
# plt.show()



##____________________________________________


# plt.scatter(length,theta,s=10)
# plt.legend(title='{}'.format(year),title_fontsize=10)
# plt.xlabel('Длина импульса, нс')
# plt.ylabel('$\Theta$,º')
# plt.show()

# # plt.scatter(length,d0,s=10)
# # plt.legend(title='{}'.format(year),title_fontsize=10)
# # plt.ylabel('d0, м')
# # plt.xlabel('Длина импульса, нс')
# # plt.show()

# plt.scatter(length,nmbr_pmt,s=10)
# plt.legend(title='{}'.format(year),title_fontsize=10)
# plt.xlabel('Длина импульса, нс')
# plt.ylabel('Число ФЭУ')
# plt.show()

# plt.scatter(ratio_or_max,theta,s=10)
# plt.legend(title='{}'.format(year),title_fontsize=10)
# plt.ylabel('Ratio')
# plt.xlabel('Число ФЭУ')
# # plt.xlim(0,50)
# plt.show()

# plt.scatter(length,a2,s=10)
# plt.legend(title='{}'.format(year),title_fontsize=10)
# plt.xlabel('d0, м')
# plt.ylabel('a2')
# # plt.xlim(0,50)
# plt.show()


# # ##____________________________________________

# # ###  PARAMETRS OF THE FRONT 

# color='blue'
# x_a0_a1_a2=np.array(x_a0)-4.814323473013257

# plt.scatter(x_a0_a1_a2,y_a0,color=color)
# plt.errorbar(x_a0_a1_a2, y_a0, yerr = y_a0_std,fmt ='o',markersize=0.1, capsize=3,color=color)
# plt.ylabel('${a_0}$',fontsize=22)
# plt.xlabel('d, м',fontsize=22)
# plt.ylim(0,max(y_a0)+250)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# # plt.xlim(None,220)

# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_a0_d0.pdf',bbox_inches='tight')
# plt.show()

# # # x_a1_elec,y_a1_elec,y_a1_std_elec=x_a1,y_a1,y_a1_std

# # # plt.scatter(np.array(x_a1_elec)+val,y_a1_elec,color='darkorange')
# # # plt.errorbar(np.array(x_a1_elec)+val, y_a1_elec, yerr = y_a1_std_elec,fmt ='o',markersize=0.1, capsize=3,color='darkorange')
    
# plt.scatter(x_a0_a1_a2,y_a1,color=color)
# plt.errorbar(x_a0_a1_a2, y_a1, yerr = y_a1_std,fmt ='o',markersize=0.1, capsize=3,color=color)
# plt.ylabel('${a_1}$',fontsize=22)
# plt.xlabel('d, м',fontsize=22)
# # plt.ylim(0,5)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.ylim(min(y_a1)-0.7,max(y_a1)+0.7)

# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_a1_d0.pdf',bbox_inches='tight')
# plt.show()


# # x_a2_elec,y_a2_elec,y_a2_std_elec=x_a2,y_a2,y_a2_std

# # plt.scatter(np.array(x_a2_elec),y_a2_elec,color=color)
# # plt.errorbar(np.array(x_a2_elec), y_a2_elec, yerr = y_a2_std_elec,fmt ='o',markersize=0.1, capsize=3,color=color)

# plt.scatter(x_a0_a1_a2,y_a2,color=color)
# plt.errorbar(x_a0_a1_a2, y_a2, yerr = y_a2_std,fmt ='o',markersize=0.1, capsize=3,color=color)
# plt.ylabel('${a_2}$',fontsize=22)
# plt.xlabel('d, м',fontsize=22)
# # plt.xlim(None,220)
# plt.ylim(0,max(y_a2)+0.003)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_a2_d0.pdf',bbox_inches='tight')
# plt.show()

# #___________________________

# # DEPENDANCES

# #theta (length)

# plt.scatter(x_th_length,np.array(y_th_length)*180/np.pi,color='dodgerblue')
# plt.errorbar(x_th_length, np.array(y_th_length)*180/np.pi, yerr = np.array(y_th_length_std)*180/np.pi,fmt ='o',markersize=0.1, capsize=3,color='dodgerblue')
# plt.ylabel('$\Theta$, °',fontsize=16)
# plt.xlabel('Длина импульса, нс',fontsize=16)
# # plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_theta_lenght.pdf',bbox_inches='tight')

# plt.show()

# # length (nbr pmt)

# plt.scatter(x_length_pmt,y_length_pmt,color='dodgerblue')
# plt.errorbar(x_length_pmt, y_length_pmt, yerr = y_length_pmt_std,fmt ='o',markersize=0.1, capsize=3,color='dodgerblue')
# plt.xlabel('Длина импульса, нс',fontsize=16)
# plt.ylabel('Число ФЭУ',fontsize=16)
# # plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_lenght_pmt.pdf',bbox_inches='tight')

# plt.show()

# # length (nbr pmt)

# plt.scatter(x_theta_ratio,np.array(y_theta_ratio)*180/np.pi,color='dodgerblue')
# plt.errorbar(x_theta_ratio, np.array(y_theta_ratio)*180/np.pi, yerr = np.array(y_theta_ratio_std)*180/np.pi,fmt ='o',markersize=0.1, capsize=3,color='dodgerblue')
# plt.xlabel('Число ФЭУ',fontsize=22) #number of PMT, pas = 10, results[4] 
# # plt.xlabel('Длина импульса, нс',fontsize=22) #length, pas = 100, results[5] 
# # plt.xlabel('max($A_∑$) / $\overline{BG_∑}$',fontsize=22) #ratio, pas = 2, results[8] 
# # plt.xlabel('$a_0$',fontsize=22) #a0, pas = 25, results[9] 
# plt.ylabel('$\Theta$,°',fontsize=22)
# plt.ylim(0,50)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_theta_pmt.pdf',bbox_inches='tight')
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_theta_lenght.pdf',bbox_inches='tight')
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_theta_a0.pdf',bbox_inches='tight')
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/exp_theta_ratio.pdf',bbox_inches='tight')

# plt.show()



