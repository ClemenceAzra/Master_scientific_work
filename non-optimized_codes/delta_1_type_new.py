import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###=============================================================
###=============================================================

## !!! CHANGE ONLY HERE

# nuclei='Fe'
# nuclei='N'
nuclei='P'

# Choose the altitude
H=500 #m, 900 

# Choose the energy
En=30 #PeV, 10

integral=0.5

# Choose the analyse
analyse='only_sig' #elec, bg
# analyse='electronic' #elec, bg
telescope='sph2'

angle_th='axis_pmt'

if telescope=='sph3':
    nucl_or_no='/{}'.format(nuclei) #/{}/
else:
    nucl_or_no=''.format(nuclei) #/{}/

## !!! END CHANGE ONLY HERE

###=============================================================
###=============================================================

##=== Open files

#Real angles
# if telescope=='sph2':

real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT
real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) #number of PMT

#Obtain results

results=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/{}/angles/{}{}/{}_total_{}m_{}PeV_{}_{}_{}.csv'.format(telescope,analyse,nucl_or_no,analyse,H,En,nuclei,telescope,angle_th),header=None,sep='\s+')
# results=resultss[resultss[4]>60]
files=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/{}/files/{}{}/{}_files_{}m_{}PeV_{}_{}_{}.csv'.format(telescope,analyse,nucl_or_no,analyse,H,En,nuclei,telescope,angle_th),header=None)[0])
# files=list(np.array(filess)[resultss[4]>60])

theta=results[0]
phi=results[1]
 
nmbr_pmt=results[4]
lenght=results[5]
num_files=results[6]
d0=results[7]/1.6292064424266897  #0.9051146902370498 / 1.6292064424266897 
ratio_or_max=results[8]
a0,a1,a2=results[9],results[10],results[11]

# print(nuclei,np.around(np.mean(a0),2),'+-',np.around(np.std(a0),2),np.around(np.mean(a1),3),'+-',np.around(np.std(a1),4),np.around(np.mean(a2),5),'+-',np.around(np.std(a2),5))

# #SPH 3
if telescope=='sph3':

    # results=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/sph3/angles/only_sig/P/only_sig_total_1000m_10PeV_P_sph3_new.csv'.format(nuclei,H,En,nuclei,angle_th),header=None,sep='\s+')
    # files=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/files/sph3/{}/{}/only_sig_files_{}m_{}PeV_{}_sph3_{}.csv'.format(nuclei,H,En,nuclei,angle_th),header=None)[0])

    
    # theta=results[0]
    # phi=results[1]
 
    # nmbr_pmt=results[4]
    # lenght=results[5]
    # num_files=results[6]
    # d0=results[7]
    # ratio_or_max=results[8]
    # a0,a1,a2=results[9],results[10],results[11]
   
    # theta_real=15/180*np.pi
    real_phi_sph3_Fe=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_{}'.format(nuclei),header=None,sep='\s+')[0])
    theta_real=15*np.pi/180


###=============================================================
###=============================================================

##=== Calculation of errors

#Real angles for each obtain files
if telescope=='sph2':  
    # theta_real=[real_theta[int(files[i][-8:-5])-1] for i in range(0,len(files))]
    phi_real=[real_phi[int(files[i][-8:-5])-1] for i in range(0,len(files))]

    if angle_th=='05':
        theta_real=[5/180*np.pi]*len(files)
    elif angle_th=='30':
        theta_real=[30/180*np.pi]*len(files)
    else:
        theta_real=[real_theta[int(files[i][-8:-5])-1] for i in range(0,len(files))]

if telescope=='sph3':  
    phi_real=[real_phi_sph3_Fe[int(files[i][-8:-5])] for i in range(0,len(files))]


#Functions

#DELTA

def f_delta(theta_fit,phi_fit,theta_real,phi_real):
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    return (np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit))*180/np.pi

#DEPENDANCES

#Dependance delta ()

def analyse_dep(results,colonne,step):
    table=results[[colonne,'delta_max']] #dataframe of nbr PMT (col 1) and delta (col 2)
    lenght=np.arange(results[colonne].min(),results[colonne].max()+5,step)
    x_step=[]
    y_mean=[]
    y_std=[]
    for i in range(0,len(lenght)-1):
        y=table[ (table[colonne]>lenght[i]) & (table[colonne]<lenght[i+1]) ]['delta_max']
        if len(y)>2:
            x_step.append(lenght[i])
            y_mean.append(y.mean())
            y_std.append(y.std())
    return  x_step, y_mean, y_std

#Dependance front param (d0)

def analyse_param(results,colonne,step):
    colonne_d0=7
    table=results[[colonne,colonne_d0]] #dataframe of param a (col 1) and d0 (col 2)
    lenght=np.arange(20,results[colonne_d0].max()+5,step)
    x_step=[]
    y_mean=[]
    y_std=[]
    for i in range(0,len(lenght)-1):
        y=table[ (table[colonne_d0]>lenght[i]) & (table[colonne_d0]<lenght[i+1]) ][colonne]
        if len(y)>2:
            x_step.append(lenght[i])
            y_mean.append(y.mean())
            y_std.append(y.std())
    return  x_step, y_mean, y_std
#==========

delta_max=f_delta(theta,phi,theta_real,phi_real)

print(np.mean(delta_max),'+-',np.std(delta_max))
print(nuclei, H)

#=================================
#=================================


results['delta_max']=delta_max

results['real_th']=np.array(theta_real)*180/np.pi

results_limits=results[ (results[5]>100) & (results['delta_max']>15)]

bad_files=np.array(files)[np.where(delta_max<1)]
perc_below_1=len(bad_files)*100/len(files)

bad_files=delta_max[(delta_max<3) & (delta_max>1)]
perc_below_3=len(bad_files)*100/len(files)

bad_files=delta_max[(delta_max<5) & (delta_max>3)]
perc_below_5=len(bad_files)*100/len(files)

bad_files=np.array(files)[np.where(delta_max>5)]
perc_below_10=len(bad_files)*100/len(files)

d0_ind_100=np.array(files)[np.where(d0<50)]


print(perc_below_1,',',perc_below_3,',',perc_below_5,',',perc_below_10)

total_files=len(files)*100/6000

print(len(results)*100/6000)

# np.savetxt('/Users/clemence/Documents/{}m_{}_{}PeV_{}.csv'.format(H,nuclei,En,angle_th),bad_files,fmt='%s')

###=============================================================
###=============================================================

## RESULTS OF FUNCTION

#Dependance: delta(nbr PMT)

delta_pmt=analyse_dep(results,4,2)
x_delta_pmt,y_delta_pmt,y_delta_pmt_std=delta_pmt[0],delta_pmt[1],delta_pmt[2]

#Dependance: delta(lenght of the signal)

delta_lenght=analyse_dep(results,5,50)
x_delta_lenght,y_delta_lenght,y_delta_lenght_std=delta_lenght[0],delta_lenght[1],delta_lenght[2]

#Dependance: delta(d0)

delta_d0=analyse_dep(results,7,20)
x_delta_d0,y_delta_d0,y_delta_d0_std=delta_d0[0],delta_d0[1],delta_d0[2]


#Dependance: delta(ratio or max)
delta_ratio_or_max=analyse_dep(results,8,2)
x_delta_ratio,y_delta_ratio,y_delta_ratio_std=delta_ratio_or_max[0],delta_ratio_or_max[1],delta_ratio_or_max[2]


#Dependance:(real theta)

delta_theta=analyse_dep(results,'real_th',1)
x_delta_theta,y_delta_theta,y_delta_theta_std=delta_theta[0],delta_theta[1],delta_theta[2]

#________

#Dependance:(a0)

delta_a0=analyse_param(results,9,20)
x_a0,y_a0,y_a0_std=delta_a0[0],delta_a0[1],delta_a0[2]

#Dependance:(a1)

delta_a1=analyse_param(results,10,20)
x_a1,y_a1,y_a1_std=delta_a1[0],delta_a1[1],delta_a1[2]

#Dependance:(a2)

delta_a2=analyse_param(results,11,20)
x_a2,y_a2,y_a2_std=delta_a2[0],delta_a2[1],delta_a2[2]



###=============================================================
###=============================================================

##==== GRAPHICS

params = {"xtick.direction": "in", "ytick.direction": "in"} #ticks in box
plt.rcParams.update(params)

####===== 


# plt.hist(delta_max,label='{} \n< $\delta$ > ={:.2f} \n$\delta$ max={:.2f}º'.format(nuclei,np.mean(delta_max),max(delta_max)))
# plt.xlabel('$\delta$')
# plt.legend(title='H={} m, E={} PeV  \nreal th={}'.format(H,En,angle_th))
# # plt.yscale('log')
# plt.show()

plt.scatter(nmbr_pmt,delta_max,s=1)
plt.ylabel('$\delta$')
plt.xlabel('Число ФЭУ')
plt.show()


plt.scatter(d0,delta_max,s=1)
plt.ylabel('$\delta$, º')
plt.xlabel('d0, м')
plt.show()

plt.scatter(ratio_or_max,delta_max,s=1)
plt.ylabel('$\delta$, º')
plt.xlabel('ratio')
plt.show()

# plt.scatter(lenght,delta_max,s=1)
# plt.ylabel('$\delta$, º')
# plt.xlabel('lenght, ns')
# plt.show()

# plt.scatter(d0,a0,s=1)
# plt.xlabel('d0, м')
# # plt.ylim(-0.7,-0.0)
# plt.ylim(10,150)
# plt.legend(title='{} - a0'.format(analyse),title_fontsize=20)
# plt.show()

###=============================================================

## DEPENDANCES

#Dependance: delta(nbr PMT)

# plt.scatter(x_delta_pmt,y_delta_pmt,color='grey')
# plt.errorbar(x_delta_pmt, y_delta_pmt, yerr = y_delta_pmt_std,fmt ='o',markersize=0.1, capsize=3,color='grey')
# plt.ylabel('<$\delta$>, º',fontsize=14)
# plt.xlabel('Число ФЭУ',fontsize=14)
# plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()

# # Dependance: delta(lenght of the signal)

# plt.scatter(x_delta_lenght,y_delta_lenght,color='dodgerblue')
# plt.errorbar(x_delta_lenght, y_delta_lenght, yerr = y_delta_lenght_std,fmt ='o',markersize=0.1, capsize=3,color='dodgerblue')
# plt.ylabel('<$\delta$>, º',fontsize=14)
# plt.xlabel('Длина импульса, нс ',fontsize=14)
# plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()


# plt.scatter(x_delta_ratio,y_delta_ratio,color='green')
# plt.errorbar(x_delta_ratio, y_delta_ratio, yerr = y_delta_ratio_std,fmt ='o',markersize=0.1, capsize=3,color='green')
# plt.ylabel('<$\delta$>, º',fontsize=14)
# plt.xlabel('Отношение макс. сигнал / фон ',fontsize=14)
# plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()

# #Dependance: delta(d0)

# plt.scatter(x_delta_d0,y_delta_d0,color='orange')
# plt.errorbar(x_delta_d0, y_delta_d0, yerr = y_delta_d0_std,fmt ='o',markersize=0.1, capsize=3,color='orange')
# plt.ylabel('<$\delta$>, º',fontsize=14)
# plt.xlabel('Ось ливня, м',fontsize=14)
# plt.ylim(0,5)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()


# #Dependance: delta(theta)
# if analyse!='sph3':
#     plt.scatter(x_delta_theta,y_delta_theta,color='purple')
#     plt.errorbar(x_delta_theta, y_delta_theta, yerr = y_delta_theta_std,fmt ='o',markersize=0.1, capsize=3,color='purple')
#     plt.ylabel('<$\delta$>, º',fontsize=14)
#     plt.xlabel('$\Theta$r, º',fontsize=14)
#     plt.ylim(0,5)
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     plt.show()

# #______

if analyse=='electronic':
    color='orange'
else:
    color='mediumvioletred'

val=3
color='crimson'

# x_a0_elec,y_a0_elec,y_a0_std_elec=x_a0,y_a0,y_a0_std

# plt.scatter(np.array(x_a0_elec)+val,y_a0_elec,color='darkorange')
# plt.errorbar(np.array(x_a0_elec)+val, y_a0_elec, yerr = y_a0_std_elec,fmt ='o',markersize=0.1, capsize=3,color='darkorange')
    
plt.scatter(x_a0,y_a0,color=color)
plt.errorbar(x_a0, y_a0, yerr = y_a0_std,fmt ='o',markersize=0.1, capsize=3,color=color)
plt.ylabel('${a_0}$',fontsize=22)
plt.xlabel('d, м',fontsize=22)
# plt.ylim(0,5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/a0_d0.pdf',bbox_inches='tight')
plt.show()

# x_a1_elec,y_a1_elec,y_a1_std_elec=x_a1,y_a1,y_a1_std

# plt.scatter(np.array(x_a1_elec)+val,y_a1_elec,color='darkorange')
# plt.errorbar(np.array(x_a1_elec)+val, y_a1_elec, yerr = y_a1_std_elec,fmt ='o',markersize=0.1, capsize=3,color='darkorange')
    
plt.scatter(x_a1,y_a1,color=color)
plt.errorbar(x_a1, y_a1, yerr = y_a1_std,fmt ='o',markersize=0.1, capsize=3,color=color)
plt.ylabel('${a_1}$',fontsize=22)
plt.xlabel('d, м',fontsize=22)
# plt.ylim(0,5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/a1_d0.pdf',bbox_inches='tight')
plt.show()


# x_a2_elec,y_a2_elec,y_a2_std_elec=x_a2,y_a2,y_a2_std

# plt.scatter(np.array(x_a2_elec)+val,y_a2_elec,color='darkorange')
# plt.errorbar(np.array(x_a2_elec)+val, y_a2_elec, yerr = y_a2_std_elec,fmt ='o',markersize=0.1, capsize=3,color='darkorange')

plt.scatter(x_a2,y_a2,color=color)
plt.errorbar(x_a2, y_a2, yerr = y_a2_std,fmt ='o',markersize=0.1, capsize=3,color=color)
plt.ylabel('${a_2}$',fontsize=22)
plt.xlabel('d, м',fontsize=22)
# plt.ylim(0,5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/a2_d0.pdf',bbox_inches='tight')

plt.show()

#==========================================================

#PARAMETERS OF THE FRONT AS A FUNCTION OF Z

# z=[1,7,26] 
# parama0=[60.98016394146646,60.595925558827645,61.01309282862487]
# parama1=[-0.28632265800890194,-0.28802335970382986,-0.28855796178754106]
# parama2=[0.0060999650633440095,0.006090823574075287 ,0.00607827531170087]

# plt.scatter(z,parama0,color='blue')
# plt.errorbar(z, parama0, yerr = [21.411200542951615,21.01765992383334,21.018458976255644],fmt ='o',markersize=0.1, capsize=3,color='blue')
# plt.xlabel('Z',fontsize=16)
# plt.ylabel('<a0>',fontsize=16)
# plt.show()

# plt.scatter(z,parama1,color='darkred')
# plt.errorbar(z, parama1, yerr = [0.09201208041945787,0.08646478818625497,0.08389814384086221],fmt ='o',markersize=0.1, capsize=3,color='darkred')
# plt.xlabel('Z',fontsize=16)
# plt.ylabel('<a1>',fontsize=16)
# plt.show()

# plt.scatter(z,parama2,color='orange')
# plt.errorbar(z, parama2, yerr = [0.00042755649386061874,0.0004152277031658902,0.0004109944842214194],fmt ='o',markersize=0.1, capsize=3,color='orange')
# plt.xlabel('Z',fontsize=16)
# plt.ylabel('<a2>',fontsize=16)
# plt.show()




