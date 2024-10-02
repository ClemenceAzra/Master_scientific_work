import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

axis_results = pd.read_csv('/Users/clemence/Documents/Магистратура_наука/Научная_работа/data_axis.csv')
delta = pd.read_csv('/Users/clemence/Documents/Магистратура_наука/Научная_работа/delta.csv')

params = {"xtick.direction": "in", "ytick.direction": "in", 'axes.axisbelow':True} #ticks in box
plt.rcParams.update(params)

fig, axs= plt.subplots(1, 2, sharex='col', gridspec_kw={'hspace': 0, 'wspace': 0.2})
fig.set_figwidth(10) 
fig.set_figheight(3.5) 
(ax1), (ax2) = axs

width_bins=np.arange(0,4.5,0.25)
ax1.hist(delta['delta_1000m'], bins=width_bins,color='orange',histtype='step',label='1000 m', linewidth=1.3)
ax1.hist(delta['delta_500m'], bins=width_bins,color='green',histtype='step',label='500 m', linewidth=1.3)
ax1.set_xlabel('$\Omega$, deg', fontsize = 18)
ax1.set_ylabel('N$_{events}$', fontsize = 18)
ax1.set_xlim(0, 4.5)
ax1.legend(fontsize=15)
ax1.set_ylim(0, 1300)
ax1.tick_params(labelsize = 15)
ax1.text(0.15, 1150, '$\it{(a)}$', fontsize = 20)

nuclei = 'p'
ax2.scatter(axis_results['1000m_x'],axis_results['1000m_mean'],color='orange',label='1000 m', edgecolor ='black',s=15)
ax2.errorbar(axis_results['1000m_x'],axis_results['1000m_mean'], yerr = axis_results['1000m_std'],fmt ='o',markersize=1, capsize=3,color='orange')
ax2.scatter(axis_results['500m_x'],axis_results['500m_mean'],color='green',label='500 m', edgecolor ='black',s=15)
ax2.errorbar(axis_results['500m_x'],axis_results['500m_mean'], yerr = axis_results['500m_std'],fmt ='o',markersize=1, capsize=3,color='green')
ax2.set_ylabel('∆d, m',fontsize=18)
ax2.set_xlabel('d, m',fontsize=18)
ax2.legend(loc='upper center',fontsize=15,title_fontsize=15)
ax2.set_ylim(0, 40)
ax2.tick_params(labelsize = 15)
ax2.text(-5., 35, '$\it{(b)}$', fontsize = 20)

plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/ЯФ-2024/axis_delta_sph3.pdf',bbox_inches='tight')

plt.show()


