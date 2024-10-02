import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

read_file = pd.read_csv('/Users/clemence/Documents/Магистратура_наука/all_dbg_params_for_events.txt', sep = '\s+')
read_file_2013 = pd.read_csv('/Users/clemence/Documents/Магистратура_наука/2013_events_with_telemetry.csv')

read_file_angles_algo = pd.read_csv('/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/RESULTS/ALL_RESULTS/experiment/experiment_files_both_years_2_100.csv', names = ['files'])

#_______

#Inclination

def inclination(enter_file, year, clin1, clin2, name_num_ev, val_plus, val_minus):

    data_2012_2013 = enter_file[enter_file['year'].isin([int(year)])]
        
    rd1 = np.array(data_2012_2013[clin1]+val_plus)
    rd2 = np.array(data_2012_2013[clin2]-val_minus)
    
    
       # // calc_mosaic_incline_vector(rd1, rd2)
    rd1 *= np.pi / 180
    rd2 *= np.pi / 180
    
    xa = np.tan(rd1)
    xy = np.tan(rd2)
    
    rdres = 1.0 + xa*xa + xy*xy
    rdres = np.sqrt(rdres)
    rdres = 1.0 / rdres
    rdres = np.arccos(rdres)
    rdres *= 180./np.pi
    
    data_event = pd.DataFrame({'event': list(data_2012_2013[name_num_ev]), 'inclination': rdres, 'year':list(data_2012_2013['year'])})
    data_event['inclination'] = list(data_event['inclination'].ffill())
    
    #_______
    
    #Angles of analyse
    
    files_name = list(map(int, read_file_angles_algo['files'].str.slice(11, 16)))
    
    inclination_values = data_event[data_event['event'].isin(files_name)].dropna().groupby('event').mean()

    return inclination_values

inclin_2012 = inclination(read_file, '2012', 'Clin1', 'Clin2', 'kevent', 4.4, 2.6)
inclin_2013 = inclination(read_file_2013, '2013', 'Clin1_y', 'Clin2_y', 'Eid', 0, 0)

plt.hist([list(inclin_2012['inclination']), list(inclin_2013['inclination'])], 
          bins = np.arange(0, 8), edgecolor = 'black', stacked = True, 
          color = ['orange', 'blue'], label = ['2012', '2013'])
plt.xlabel('$\Delta$, º', fontsize = 16)
plt.ylabel('Количество событий', fontsize = 16)
plt.grid(True)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.xlim(0)
plt.legend()
# plt.savefig('/Users/clemence/Documents/Магистратура_наука/Мои статьи/Статья_3_УЗ/Figs/inclination.eps', bbox_inches = 'tight')
plt.show()


