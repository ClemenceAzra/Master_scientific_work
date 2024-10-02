import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

telescope='SPHERE-3'

###=============================================================
###=============================================================
dirname='/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/'

path_only_sig=dirname+'RESULTS/ALL_RESULTS/sph2/angles/only_sig'
path_sig_bg=dirname+'RESULTS/ALL_RESULTS/sph2/angles/bg_sig'
path_elec=dirname+'RESULTS/ALL_RESULTS/sph2/angles/electronic'

path_sph3_Fe=dirname+'RESULTS/ALL_RESULTS/sph3/angles/only_sig/Fe'
path_sph3_N=dirname+'RESULTS/ALL_RESULTS/sph3/angles/only_sig/N'
path_sph3_P=dirname+'RESULTS/ALL_RESULTS/sph3/angles/only_sig/P'


real_theta=np.array(pd.read_csv(dirname+'Initial_data/angles_theta',header=None,sep='\s+')[0]) #number of PMT
real_phi=np.array(pd.read_csv(dirname+'Initial_data/angles_phi',header=None,sep='\s+')[0]) #number of PMT

real_phi_sph3_Fe=np.array(pd.read_csv(dirname+'Initial_data/angels_phi_Fe_sph3',header=None,sep='\s+')[0])
real_phi_sph3_N=np.array(pd.read_csv(dirname+'Initial_data/angels_phi_N_sph3',header=None,sep='\s+')[0])
real_phi_sph3_P=np.array(pd.read_csv(dirname+'Initial_data/angels_phi_P_sph3',header=None,sep='\s+')[0])


##=== Open directories

name_files_only_sig = glob.glob(os.path.join(path_only_sig,'*'))
name_files_sig_bg = glob.glob(os.path.join(path_sig_bg,'*'))
name_files_elec = glob.glob(os.path.join(path_elec,'*'))

name_files_sph3_Fe = glob.glob(os.path.join(path_sph3_Fe,'*'))
name_files_sph3_N = glob.glob(os.path.join(path_sph3_N,'*'))
name_files_sph3_P = glob.glob(os.path.join(path_sph3_P,'*'))

###=============================================================
###=============================================================


def f_delta(theta_fit,phi_fit,theta_real,phi_real):
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    delta=(np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit))*180/np.pi
    return  delta, theta_fit, theta_real

deltas_bar=np.arange(0,20+0.5,0.5)

def graph(delta):
    delta_perc=[]
    for i in range(0,len(deltas_bar)-1):
        # delta_df=pd.DataFrame({'delta':delta_max})
        delta_perc.append(len(delta[ (delta>deltas_bar[i]) & (delta<deltas_bar[i+1])])*100/len(delta))
    return delta_perc


def results(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    nbr_pmt=[]
    length=[]
    
    theta_real=[]
    theta_fit=[]
    phi_fit=[]
    phi_real=[]
    
    for files in range(0,len(name_files)):
        # print(name_files[files])
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi[np.array(result_angle[6].astype('int')-1)]
        real_theta_analyse=real_theta[np.array(result_angle[6].astype('int')-1)]
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)[0]
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        nbr_pmt.append(result_angle[4])
        d0.append(result_angle[7])
        ratio.append(result_angle[8])
        length.append(result_angle[5])
        
        theta_real.append(real_theta_analyse)
        phi_real.append(real_phi_analyse)
        theta_fit.append(list(result_angle[0]))
        phi_fit.append(list(result_angle[1]))

        
    return delta_total,delta_multiple,d0,ratio,nbr_pmt,length , theta_real, theta_fit, phi_real, phi_fit
 
def results_sph2(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    nbr_pmt=[]
    length=[]
    theta_real=[]
    theta_fit=[]
    phi_fit=[]
    phi_real=[]
    for files in range(0,len(name_files)):
        # print(name_files[files])
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi[np.array(result_angle[6].astype('int')-1)]
        real_theta_analyse=real_theta[np.array(result_angle[6].astype('int')-1)]
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)[0]
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        nbr_pmt.append(result_angle[4])
        d0.append(result_angle[7])
        ratio.append(result_angle[8])
        length.append(result_angle[5])
         
        theta_real.append(real_theta_analyse)
        theta_fit.append(list(result_angle[0]))
        phi_real.append(real_phi_analyse)
        phi_fit.append(list(result_angle[1]))

                
    return delta_total,delta_multiple,d0,ratio,nbr_pmt,length, theta_real, theta_fit, phi_real, phi_fit
    
def ratio_all(delta,ratio):
    delta_ratio=pd.DataFrame({'delta':delta,'ratio':ratio})
    delta_ratio['ratio']=delta_ratio['ratio']. round(0)

    moy_delta=delta_ratio.groupby(delta_ratio['ratio']).mean()
    std_delta=delta_ratio.groupby(delta_ratio['ratio']).std()
    
    x=list(moy_delta.index)
    mean=list(moy_delta['delta'])
    std=list(std_delta['delta'])
    
    return x,mean,std

#=====================================================================================
#=====================================================================================

def results_sph3_Fe(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    nbr_pmt=[]
    length=[]
    theta_real=[]
    theta_fit=[]
    phi_fit=[]
    phi_real=[]
    for files in range(0,len(name_files)):
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi_sph3_Fe[np.array(result_angle[6].astype('int'))]
        real_theta_analyse=[15/180*np.pi]*len(real_phi_analyse)
        
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)[0]
        
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        nbr_pmt.append(result_angle[4])
        d0.append(result_angle[7])
        ratio.append(result_angle[8])
        length.append(result_angle[5])
        theta_real.append(real_theta_analyse)
        theta_fit.append(list(result_angle[0]))
        phi_fit.append(list(result_angle[1]))
        phi_real.append(real_phi_analyse)

        


        
    return delta_total,delta_multiple,d0,ratio,nbr_pmt,length,theta_real, theta_fit, phi_real, phi_fit

def results_sph3_N(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    nbr_pmt=[]
    length=[]
    theta_real=[]
    theta_fit=[]
    phi_fit=[]
    phi_real=[]
    for files in range(0,len(name_files)):
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi_sph3_N[np.array(result_angle[6].astype('int'))]
        real_theta_analyse=[15/180*np.pi]*len(real_phi_analyse)
        
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)[0]
        
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        nbr_pmt.append(result_angle[4])
        d0.append(result_angle[7])
        ratio.append(result_angle[8])
        length.append(result_angle[5])
        theta_real.append(real_theta_analyse)
        theta_fit.append(list(result_angle[0]))
        phi_fit.append(list(result_angle[1]))        
        phi_real.append(real_phi_analyse)


        
    return delta_total,delta_multiple,d0,ratio,nbr_pmt,length,theta_real, theta_fit, phi_real, phi_fit

def results_sph3_P(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    nbr_pmt=[]
    length=[]
    theta_fit=[]
    phi_fit=[]
    phi_real=[]
    theta_real=[]
    for files in range(0,len(name_files)):
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi_sph3_P[np.array(result_angle[6].astype('int'))]
        real_theta_analyse=[15/180*np.pi]*len(real_phi_analyse)
        
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)[0]
        
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        nbr_pmt.append(result_angle[4])
        d0.append(result_angle[7])
        ratio.append(result_angle[8])
        length.append(result_angle[5])
        theta_real.append(real_theta_analyse)
        theta_fit.append(list(result_angle[0]))
        phi_fit.append(list(result_angle[1]))
        phi_real.append(real_phi_analyse)


        
    return delta_total,delta_multiple,d0,ratio,nbr_pmt,length,theta_real, theta_fit, phi_real, phi_fit


#Dependance delta ()

def analyse_dep(delta,col_x,step):
    table=pd.DataFrame({'delta':delta,'col_x':col_x}) #dataframe of nbr PMT (col 1) and delta (col 2)
    lenght=np.arange(table['col_x'].min(),table['col_x'].max(),step)
    x_step=[]
    y_mean=[]
    y_std=[]
    for i in range(0,len(lenght)-1):
        y=table[ (table['col_x']>lenght[i]) & (table['col_x']<lenght[i+1]) ]['delta']
        if len(y)>2:
            x_step.append(lenght[i])
            y_mean.append(y.mean())
            y_std.append(y.std())
    return  x_step, y_mean, y_std



###=============================================================
###=============================================================

#ONLY SIG

all_results_only_sig=results(name_files_only_sig)
delta_total,delta_multiple,d0,ratio,nbr_pmt,length_only_sig=all_results_only_sig[0],all_results_only_sig[1],all_results_only_sig[2],all_results_only_sig[3],all_results_only_sig[4],all_results_only_sig[5]

theta_real_total_only_sig=np.concatenate(all_results_only_sig[6])
theta_fit_total_only_sig=np.concatenate(all_results_only_sig[7])

#SIG AND BG

all_results_sig_bg=results(name_files_sig_bg)
delta_total_sig_bg,delta_multiple_sig_bg,d0_bg,ratio_bg=all_results_sig_bg[0],all_results_sig_bg[1],all_results_sig_bg[2],all_results_sig_bg[3]

#ELECTRONIC

all_results_elec=results_sph2(name_files_elec)
delta_total_elec,delta_multiple_elec,d0_elec,ratio_elec,nbr_pmt_elec,length_elec=all_results_elec[0],all_results_elec[1],all_results_elec[2],all_results_elec[3],all_results_elec[4],all_results_elec[5]

theta_real_total_elec=np.concatenate(all_results_elec[6])
theta_fit_total_elec=np.concatenate(all_results_elec[7])
# delta_total_elec=np.concatenate(delta_total_elec)

#SPH3

all_results_sph3_N=results_sph3_N(name_files_sph3_N)
delta_total_sph3_N,delta_multiple_sph3_N,d0_sph3_N,ratio_sph3_N,nbr_pmt_sph3_N,length_sph3_N,real_theta_sph3_N=all_results_sph3_N[0],all_results_sph3_N[1],all_results_sph3_N[2],all_results_sph3_N[3],all_results_sph3_N[4],all_results_sph3_N[5],all_results_sph3_N[6]

all_results_sph3_Fe=results_sph3_Fe(name_files_sph3_Fe)
delta_total_sph3_Fe,delta_multiple_sph3_Fe,d0_sph3_Fe,ratio_sph3_Fe,nbr_pmt_sph3_Fe,length_sph3_Fe,real_theta_sph3_Fe=all_results_sph3_Fe[0],all_results_sph3_Fe[1],all_results_sph3_Fe[2],all_results_sph3_Fe[3],all_results_sph3_Fe[4],all_results_sph3_Fe[5],all_results_sph3_Fe[6]

all_results_sph3_P=results_sph3_P(name_files_sph3_P)
delta_total_sph3_P,delta_multiple_sph3_P,d0_sph3_P,ratio_sph3_P,nbr_pmt_sph3_P,length_sph3_P,real_theta_sph3_P=all_results_sph3_P[0],all_results_sph3_P[1],all_results_sph3_P[2],all_results_sph3_P[3],all_results_sph3_P[4],all_results_sph3_P[5],all_results_sph3_P[6]

# deltaaa=delta_total_sph3_Fe+delta_total_sph3_Fe+delta_total_sph3_P
# sds=[np.std(deltaaa[i]) for i in range(0,len(deltaaa))]
# aa=np.mean(sds)

deltaaa=delta_total_sph3_Fe+delta_total_sph3_N+delta_total_sph3_P
theta_real_sph3=real_theta_sph3_N+real_theta_sph3_Fe+real_theta_sph3_P
ds_delta_sph3=np.mean([np.std(deltaaa[i]) for i in range(0,len(deltaaa))])
ds_sph3=np.std([len(deltaaa[i]) for i in range(0,len(deltaaa))])*100/6000

# ds=np.std([len(delta_total_elec[i]) for i in range(0,len(delta_total_elec))])*100/6000

ds_delta_only_sig=np.mean([np.mean(delta_total[i]) for i in range(0,len(delta_total))])
ds_delta_files=np.mean([len(delta_total[i]) for i in range(0,len(delta_total))])*100/6000


ds_delta_elec=np.mean([np.mean(delta_total_elec[i]) for i in range(0,len(delta_total_elec))])
ds_delta_file_elec=np.std([len(delta_total_elec[i]) for i in range(0,len(delta_total_elec))])*100/6000


delta_sph3_500m=list(delta_total_sph3_Fe[1])+list(delta_total_sph3_N[0])+list(delta_total_sph3_P[0])
delta_sph3_1000m=list(delta_total_sph3_Fe[0])+list(delta_total_sph3_N[1])+list(delta_total_sph3_P[1])

# results_delta_sph3=pd.DataFrame({'delta_500m':delta_sph3_500m,'delta_1000m':delta_sph3_1000m+[np.nan]*(len(delta_sph3_500m)-len(delta_sph3_1000m))})
# results_delta_sph3.to_csv('/Users/clemence/Documents/delta.csv')

###=============================================================
###=============================================================


#FOR DEPENDANCES delta ()

#Number of PMT

#======= SPHERE-2

#Only sig
if telescope=='SPHERE-2':

    concat_nbr_pmt_only_sig,concat_delta_only_sig,concat_d0_only_sig,concat_length_only_sig,concat_ratio_only_sig=np.concatenate(nbr_pmt),np.concatenate(delta_total),np.concatenate(d0),np.concatenate(length_only_sig),np.concatenate(ratio)
    pas_ratio=20
    pas_nbr_pmt=5
    color_elec='grey'
    color_val='blue'

#Elec
    # concat_nbr_pmt_only_sig,concat_delta_only_sig,concat_d0_only_sig,concat_length_only_sig,concat_ratio_only_sig=np.concatenate(nbr_pmt_elec),np.concatenate(delta_total_elec),np.concatenate(d0_elec),np.concatenate(length_elec),np.concatenate(ratio_elec)
    # pas_ratio=20
    # pas_nbr_pmt=5


#======= SPHERE-3
if telescope=='SPHERE-3':
    concat_nbr_pmt_only_sig,concat_delta_only_sig,concat_d0_only_sig,concat_length_only_sig,concat_ratio_only_sig=np.concatenate(nbr_pmt_sph3_N+nbr_pmt_sph3_Fe+nbr_pmt_sph3_P),np.concatenate(delta_total_sph3_N+delta_total_sph3_P+delta_total_sph3_Fe),np.concatenate(d0_sph3_Fe+d0_sph3_P+d0_sph3_N),np.concatenate(length_sph3_N+length_sph3_Fe+length_sph3_P),np.concatenate(ratio_sph3_N+ratio_sph3_Fe+ratio_sph3_P)
    pas_ratio=50
    pas_nbr_pmt=20
    color_val='green'




delta_pmt=analyse_dep(concat_delta_only_sig,concat_nbr_pmt_only_sig,pas_nbr_pmt)
x_delta_pmt,y_delta_pmt,y_delta_pmt_std=delta_pmt[0],delta_pmt[1],delta_pmt[2]

delta_d0=analyse_dep(concat_delta_only_sig,concat_d0_only_sig,25)
x_delta_d0,y_delta_d0,y_delta_d0_std=delta_d0[0],delta_d0[1],delta_d0[2]

delta_length=analyse_dep(concat_delta_only_sig,concat_length_only_sig,100)
x_delta_length,y_delta_length,y_delta_length_std=delta_length[0],delta_length[1],delta_length[2]

delta_ratio=analyse_dep(concat_delta_only_sig,concat_ratio_only_sig,pas_ratio)
x_delta_ratio,y_delta_ratio,y_delta_ratio_std=delta_ratio[0],delta_ratio[1],delta_ratio[2]

#other

dependances_d0_ratio=analyse_dep(concat_nbr_pmt_only_sig,concat_d0_only_sig,50)
x_dependances_d0_ratio,y_dependances_d0_ratio,y_dependances_d0_ratio_std=dependances_d0_ratio[0],dependances_d0_ratio[1],dependances_d0_ratio[2]

dependances_d0_ratio=analyse_dep(concat_nbr_pmt_only_sig,concat_length_only_sig,100)
x_dependances_d0_ratio,y_dependances_d0_ratio,y_dependances_d0_ratio_std=dependances_d0_ratio[0],dependances_d0_ratio[1],dependances_d0_ratio[2]


#Only sig
dependances_theta=analyse_dep(np.concatenate(delta_total),theta_real_total_only_sig*180/np.pi,0.5)
x_dependances_theta,y_dependances_theta,y_dependances_theta_std=dependances_theta[0],dependances_theta[1],dependances_theta[2]

#Elec
# dependances_theta=analyse_dep(np.concatenate(delta_total_elec),theta_real_total_elec*180/np.pi,0.5)
# x_dependances_theta,y_dependances_theta,y_dependances_theta_std=dependances_theta[0],dependances_theta[1],dependances_theta[2]

#SPH3
# delta_sph3=np.concatenate(delta_total_sph3_Fe).tolist()+np.concatenate(delta_total_sph3_N).tolist()+np.concatenate(delta_total_sph3_P).tolist()
# real_theta_sph3=np.concatenate(real_theta_sph3_Fe).tolist()+np.concatenate(real_theta_sph3_N).tolist()+np.concatenate(real_theta_sph3_P).tolist()

# dependances_theta=analyse_dep(delta_sph3,np.array(real_theta_sph3)*180/np.pi,0.5)
# x_dependances_theta,y_dependances_theta,y_dependances_theta_std=dependances_theta[0],dependances_theta[1],dependances_theta[2]


# dependances_d0_ratio=analyse_dep(concat_length_only_sig,concat_ratio_only_sig,10)
# x_dependances_d0_ratio,y_dependances_d0_ratio,y_dependances_d0_ratio_std=dependances_d0_ratio[0],dependances_d0_ratio[1],dependances_d0_ratio[2]


###=============================================================
###=============================================================

##==== Graphics

# params = {"xtick.direction": "in", "ytick.direction": "in"}
# plt.rcParams.update(params)

#===================

#ONLY SIG

#name_files_only_sig

pd_files=pd.Series(name_files_only_sig)


#=================================================

## DEPENDANCES

#Dependance: delta(nbr PMT)

# x_delta_pmt_sph3,y_delta_pmt_sph3,y_delta_pmt_std_sph3=x_delta_pmt,y_delta_pmt,y_delta_pmt_std
# val=1.5
# plt.scatter(np.array(x_delta_pmt_sph3)+val,y_delta_pmt_sph3,color=color_elec,label='Электронный сигнал')
# plt.errorbar(np.array(x_delta_pmt_sph3)+val, y_delta_pmt_sph3, yerr = y_delta_pmt_std_sph3,fmt ='o',markersize=0.1, capsize=3,color=color_elec)

# plt.scatter(x_delta_pmt,y_delta_pmt,color=color_val,label='Чистый сигнал')
# plt.errorbar(x_delta_pmt, y_delta_pmt, yerr = y_delta_pmt_std,fmt ='o',markersize=0.1, capsize=3,color=color_val)
# plt.ylabel('${δ}$, º',fontsize=14)
# plt.xlabel('Число ФЭУ',fontsize=14)
# # plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,7)
# plt.legend()
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/elec_delta_pmt_sph3.pdf',bbox_inches='tight')
# plt.show()

# # #Dependance: delta(d0)

# # x_delta_d0_sph3,y_delta_d0_sph3,y_delta_d0_std_sph3=x_delta_d0,y_delta_d0,y_delta_d0_std

# # val=3
# # plt.scatter(np.array(x_delta_d0_sph3)+val,y_delta_d0_sph3,color=color_elec,label='Электронный сигнал')
# # plt.errorbar(np.array(x_delta_d0_sph3)+val, y_delta_d0_sph3, yerr = y_delta_d0_std_sph3,fmt ='o',markersize=0.1, capsize=3,color=color_elec)

# plt.scatter(x_delta_d0,y_delta_d0,color=color_val,label='Чистый сигнал')
# plt.errorbar(x_delta_d0, y_delta_d0, yerr = y_delta_d0_std,fmt ='o',markersize=0.1, capsize=3,color=color_val)
# plt.ylabel('${δ}$, º',fontsize=14)
# plt.xlabel('Ось ливня, м',fontsize=14)
# # plt.ylim(0,4)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,7)
# plt.legend()
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/elec_delta_d0_sph3.pdf',bbox_inches='tight')
# plt.show()


# # # Dependance: delta(lenght of the signal)

# # x_delta_length_sph3,y_delta_length_sph3,y_delta_length_std_sph3=x_delta_length,y_delta_length,y_delta_length_std

# # val=15
# # plt.scatter(np.array(x_delta_length_sph3)+val,y_delta_length_sph3,color=color_elec,label='Электронный сигнал')
# # plt.errorbar(np.array(x_delta_length_sph3)+val, y_delta_length_sph3, yerr = y_delta_length_std_sph3,fmt ='o',markersize=0.1, capsize=3,color=color_elec)

# plt.scatter(x_delta_length,y_delta_length,color=color_val,label='Чистый сигнал')
# plt.errorbar(x_delta_length, y_delta_length, yerr = y_delta_length_std,fmt ='o',markersize=0.1, capsize=3,color=color_val)
# plt.ylabel('${δ}$, º',fontsize=14)
# plt.xlabel('Длина импульса, нс ',fontsize=14)
# # plt.ylim(0,4)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,7)
# plt.legend()
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/elec_delta_length_sph3.pdf',bbox_inches='tight')
# plt.show()


# # Dependance: delta(ratio)

# x_delta_ratio_sph3,y_delta_ratio_sph3,y_delta_ratio_std_sph3=x_delta_ratio,y_delta_ratio,y_delta_ratio_std

# val=7
# plt.scatter(np.array(x_delta_ratio_sph3)+val,y_delta_ratio_sph3,color=color_elec,label='Электронный сигнал')
# plt.errorbar(np.array(x_delta_ratio_sph3)+val, y_delta_ratio_sph3, yerr = y_delta_ratio_std_sph3,fmt ='o',markersize=0.1, capsize=3,color=color_elec)

# plt.scatter(x_delta_ratio,y_delta_ratio,color=color_val,label='Чистый сигнал')
# plt.errorbar(x_delta_ratio, y_delta_ratio, yerr = y_delta_ratio_std,fmt ='o',markersize=0.1, capsize=3,color=color_val)
# plt.ylabel('${δ}$, º',fontsize=14)
# plt.xlabel('max($A_∑$) / $\overline{BG_∑}$ ',fontsize=14)
# # plt.xlim(0,200)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,7)
# plt.legend()
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/elec_delta_ratio_sph3.pdf',bbox_inches='tight')
# plt.show()


#Others


# plt.scatter(x_dependances_d0_ratio,y_dependances_d0_ratio,color='darkred')
# plt.errorbar(x_dependances_d0_ratio, y_dependances_d0_ratio, yerr = y_dependances_d0_ratio_std,fmt ='o',markersize=0.1, capsize=3,color='darkred')
# plt.ylabel('Число ФЭУ',fontsize=14)

# plt.xlabel('Длина импульса, нс ',fontsize=14)
# # plt.xlabel('max($A_∑$) / $\overline{BG_∑}$ ',fontsize=14)
# # plt.xlabel('d0, м')
# plt.ylim(0,None)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/elec_pmt_length.pdf',bbox_inches='tight')
# plt.show()

#SPHERE-2 delta (theta)

# x_dependances_theta_elec,y_dependances_theta_elec,y_dependances_theta_elec_std=x_dependances_theta,y_dependances_theta,y_dependances_theta_std
# val=0.15
# # plt.scatter(np.array(x_dependances_theta_elec)+val,y_dependances_theta_elec,color='grey')
# # plt.errorbar(np.array(x_dependances_theta_elec)+val, y_dependances_theta_elec, yerr = y_dependances_theta_elec_std,fmt ='o',markersize=0.1, capsize=3,color='grey')

# plt.scatter(x_dependances_theta,y_dependances_theta,color='blue')
# plt.errorbar(x_dependances_theta, y_dependances_theta, yerr = y_dependances_theta_std,fmt ='o',markersize=0.1, capsize=3,color='blue')

# plt.xlabel('$\Theta_{real}, º$',fontsize=14)
# plt.ylabel('${δ}$, º',fontsize=14)

# # plt.axis('scaled')
# # plt.xlim(8,22)
# plt.ylim(0,7)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/delta_real_theta.pdf',bbox_inches='tight')
# plt.show()

#SPHERE-3


# x_dependances_theta_elec,y_dependances_theta_elec,y_dependances_theta_elec_std=x_dependances_theta,y_dependances_theta,y_dependances_theta_std
# val=0.15
# plt.scatter(np.array(x_dependances_theta_elec)+val,y_dependances_theta_elec,color='grey')
# plt.errorbar(np.array(x_dependances_theta_elec)+val, y_dependances_theta_elec, yerr = y_dependances_theta_elec_std,fmt ='o',markersize=0.1, capsize=3,color='grey')

# plt.scatter(x_dependances_theta,y_dependances_theta,color='blue')
# plt.errorbar(x_dependances_theta, y_dependances_theta, yerr = y_dependances_theta_std,fmt ='o',markersize=0.1, capsize=3,color='blue')

# plt.xlabel('$\Theta_{real}, º$',fontsize=20)
# plt.ylabel('${δ}$, º',fontsize=20)

# # plt.axis('scaled')
# # plt.xlim(8,22)
# plt.ylim(0,7)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/elec_pmt_length.pdf',bbox_inches='tight')
# plt.show()

# delta_sph3_1000m_adjust=np.array(delta_sph3_1000m)-min(delta_sph3_1000m)
# delta_sph3_500m_adjust=np.array(delta_sph3_500m)-min(delta_sph3_500m)

angles_results=pd.read_csv('/Users/clemence/Documents/Магистратура_наука/Научная_работа/delta.csv')

width_bins=np.arange(0,4.5,0.25)
plt.hist(angles_results['delta_1000m'],alpha=1, bins=width_bins,color='blue',histtype='step',label='1000 m', density = True, linestyle = '--')
plt.hist(angles_results['delta_500m'],alpha=1, bins=width_bins,color='purple',histtype='step',label='500 m', density = True)
plt.xlabel('$\Omega$, º',fontsize=18)
plt.ylabel('N',fontsize=18)
plt.xlim(0,5)
plt.legend(title='SPHERE-3', fontsize=15, title_fontsize = 15)
plt.rcParams['axes.axisbelow'] = True
plt.grid(True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(None, 0.72)
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Постер ECCR/Figures/angles_sph3.pdf',bbox_inches='tight')
plt.show()

width_bins=np.arange(0,4.5,0.25)
# plt.hist(angles_results['delta_1000m'], bins=width_bins,color='orange',histtype='step',label='1000 m, $\Omega$={:.2f}±{:.2f}'.format(np.mean(angles_results['delta_1000m']),np.std(angles_results['delta_1000m'])), linewidth=1.3)
# plt.hist(angles_results['delta_500m'], bins=width_bins,color='green',histtype='step',label='500 m, $\Omega$={:.2f}±{:.2f}'.format(np.mean(angles_results['delta_500m']),np.std(angles_results['delta_500m'])), linewidth=1.3)
# plt.xlabel('$\Omega$, º')
# plt.ylabel('Number of events')
# plt.xlim(0,4.5)
# plt.legend()
# plt.rcParams['axes.axisbelow'] = True
# plt.grid(True)
# # plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/ЯФ-2024/angles_sph3_only_sig.pdf',bbox_inches='tight')
# plt.show()

delta_500m=list(delta_total[1])+list(delta_total[2])+list(delta_total[5])+list(delta_total[7])+list(delta_total[8])+list(delta_total[9])
delta_900m=list(delta_total[0])+list(delta_total[3])+list(delta_total[4])+list(delta_total[6])+list(delta_total[10])+list(delta_total[11])
angles_sph2 = pd.DataFrame({'delta_500m':delta_500m+[np.nan]*(len(delta_900m) - len(delta_500m)), 'delta_900m': delta_900m})
# angles_sph2.to_csv('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Известия_РАН/angles_sph2.csv')

# plt.hist(delta_900m,alpha=1, bins=width_bins,color='blue',histtype='step',label='900 m, $\delta$={:.1f}±{:.1f}'.format(np.mean(delta_900m),np.std(delta_900m)), density = True)
plt.hist(delta_900m,alpha=1, bins=width_bins,color='blue',histtype='step',label='900 m', density = True, linestyle = '--')
plt.hist(delta_500m,alpha=1, bins=width_bins,color='purple',histtype='step',label='500 m', density = True)
plt.xlabel('$\Omega$, º',fontsize=18)
plt.ylabel('N',fontsize=18)
plt.xlim(0,5)
plt.ylim(None, 0.72)
plt.legend(title='SPHERE-2', fontsize=15, title_fontsize = 15)
plt.rcParams['axes.axisbelow'] = True
plt.grid(True)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Постер ECCR/Figures/angles_sph2.pdf',bbox_inches='tight')
plt.show()



# width_bins=np.arange(0,4.5,0.25)
# delta_500m_elec=list(delta_total_elec[0])+list(delta_total_elec[2])+list(delta_total_elec[3])+list(delta_total_elec[6])+list(delta_total_elec[7])+list(delta_total_elec[10])
# delta_900m_elec=list(delta_total_elec[1])+list(delta_total_elec[4])+list(delta_total_elec[5])+list(delta_total_elec[8])+list(delta_total_elec[9])+list(delta_total_elec[11])

# plt.hist(delta_900m_elec,alpha=1, bins=width_bins,color='blue',histtype='step',label='900 m, $\delta$={:.2f}±{:.2f}'.format(np.mean(delta_900m_elec),np.std(delta_900m_elec)))
# plt.hist(delta_500m_elec,alpha=1, bins=width_bins,color='purple',histtype='step',label='500 m, $\delta$={:.2f}±{:.2f}'.format(np.mean(delta_500m_elec),np.std(delta_500m_elec)))
# plt.xlabel('$\delta$, º',fontsize=18)
# plt.ylabel('Число событий',fontsize=18)
# plt.xlim(0,5)
# plt.legend(title='СФЕРА-2, «электронный сигнал»', fontsize=13, title_fontsize = 13, loc = 'upper right')
# plt.rcParams['axes.axisbelow'] = True
# plt.grid(True)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# # plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Статья_3_УЗ/Figs/angles_sph2_electronic.eps',bbox_inches='tight')
# plt.show()

# width_bins=np.arange(0,30,2.5)
# plt.hist(np.concatenate(all_results_only_sig[6]+all_results_elec[6]+all_results_sph3_N[6]+all_results_sph3_P[6]+all_results_sph3_Fe[6])*180/np.pi, bins=width_bins, alpha=0.6, color='blue', label='real', histtype='step')
# plt.hist(np.concatenate(all_results_only_sig[7]+all_results_elec[7]+all_results_sph3_N[7]+all_results_sph3_P[7]+all_results_sph3_Fe[7])*180/np.pi, bins=width_bins, alpha=0.6, color='green', label='fit', histtype='step')
# plt.legend(fontsize=16)
# plt.rcParams['axes.axisbelow'] = True
# plt.xlabel('$\Theta$, º',fontsize=18)
# plt.grid(True)
# plt.ylabel('Число событий',fontsize=18)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.xlim(0,None)
# # plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Статья_3_УЗ/Figs/all_theta_model.eps',bbox_inches='tight')
# plt.show()



# width_bins=np.arange(0,360,30)
# real_phi = np.concatenate(all_results_only_sig[8]+all_results_elec[8]+all_results_sph3_N[8]+all_results_sph3_P[8]+all_results_sph3_Fe[8])*180/np.pi
# results_phi = np.concatenate(all_results_only_sig[9]+all_results_elec[9]+all_results_sph3_N[9]+all_results_sph3_P[9]+all_results_sph3_Fe[9])*180/np.pi
# plt.hist(real_phi, bins=width_bins, alpha=0.6, color='blue', label='real', histtype='step')
# plt.hist(results_phi, bins=width_bins, alpha=0.6, color='green', label='fit', histtype='step')
# plt.legend(fontsize=16)
# plt.ylabel('Число событий',fontsize=18)
# plt.xlabel('$\phi$, º',fontsize=18)
# plt.rcParams['axes.axisbelow'] = True
# plt.grid(True)
# plt.xlim(0,None)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# # plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Статья_3_УЗ/Figs/all_phi_model.eps',bbox_inches='tight')
# plt.show()


# errors = np.histogram(real_phi, width_bins)[0] / np.histogram(results_phi, width_bins)[0]




