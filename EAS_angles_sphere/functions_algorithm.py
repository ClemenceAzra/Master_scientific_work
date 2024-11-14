import pandas as pd
import numpy as np
from iminuit import Minuit

class Functions_algorithm:
    'This class has each function necessary to determine the angle'
    def __init__(self,initial_values,files):
        self.initial_values=initial_values
        self.files=files
        self.event=self.files.event
        self.pied=self.initial_values.pied
        self.cells=self.initial_values.cells
        self.type_analyse=self.initial_values.type_analyse
        self.telescope=self.initial_values.telescope
        self.path=self.files.path
        self.altitude=self.initial_values.altitude
        self.x_mos=self.initial_values.x_mos
        self.y_mos=self.initial_values.y_mos
        self.total_number_of_pmt=self.initial_values.total_number_of_pmt
        self.margin_left  = self.initial_values.margin_left
        self.margin_right = self.initial_values.margin_right
        self.N_std_background_sum = self.initial_values.N_std_background_sum
        self.sum_of_impulse=self.sum_of_impulse()
        self.intensity_signal=self.intensity()
        self.diapason=self.diapason()
        self.front_snow=self.translation_snow_mos()
        self.amplification=self.amplification()
        self.noise_pmt=self.noise_pmt()
        self.front_up_noise=self.front_up_noise()
        self.DFS=self.DFS()
        self.neighbors=self.neighbors()
        self.angles=self.angles()
        self.lenght_impulse=self.length()
        self.axis=self.axis()

    def sum_of_impulse(self): 
        'Return the summed impulse in each channel'
        return self.event.sum(axis=1) #summarize event
    
    def intensity(self):
        'Return the intensity (max A / <BG>) for != only signal, and the maximum amplitude'   
        summed_impulse=self.sum_of_impulse #summed event
        max_signal=summed_impulse[1:-1].max() #amplitude max
        noise_pied=summed_impulse[1:self.pied][summed_impulse[1:self.pied]>0] #piedestal
        if self.type_analyse != 'only_sig':
            intensity_signal=(max_signal)/(noise_pied.mean()+ self.N_std_background_sum*noise_pied.std())
            return max_signal, intensity_signal
        return max_signal, 0
           
    def diapason(self): #diapason of the event 
       'Return the time diapason of impulse: N photons and time'
       summed_impulse=self.sum_of_impulse #summed event
       imax = summed_impulse[:-1].idxmax() #index of the max impulse --> delete the last cell -> > max of impulse
       start, stop = imax, imax
       maxpart = 1
       background=np.mean(summed_impulse[1:self.pied])+ self.N_std_background_sum*np.std(summed_impulse[1:self.pied])
       bgd = maxpart * background
       while summed_impulse[start] > bgd and start!=summed_impulse.index[0]:
           start -= 1
       while summed_impulse[stop] > bgd and stop!=summed_impulse.index[-1]:
           stop += 1
       stop  += self.margin_right
       start -= self.margin_left      
       event_diapason=self.event[start:stop] #diapason of the impulse - n photons
       cells_diapason=self.cells[start:stop] #diapason of the impulse - time       
       event_diapason.index=np.arange(len(event_diapason)) #reindex new finded diapason
       return event_diapason, cells_diapason, background
        
    def translation_snow_mos(self): 
        'Return the mosaic front (mm) on the snow (m): x, y, t'
        High =self.altitude        
        if self.type_analyse =='experiment':
                High = float(list(pd.read_csv(self.path).iloc[4])[0]) - 456  #Altitude of the telescope. 456 = Baikal altitude
        b=self.initial_values.translation_param*(High/500)  #coeff for translate x mos --> x snow
        x_snow=-b*self.x_mos                  #x on snow
        y_snow=-b*self.y_mos                  #y on snow
        t_path=((np.sqrt(High**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/self.initial_values.c_ns) #path from mosaic to snow in ns
        return x_snow, y_snow, t_path 
        
    def axis(self): 
        'Return the axis of the EAS on snow [x (m), y (m), distance from center (mm)] by 3 different methods:'
        'if the maximum is on the last circle, this is the axis'
        'if the maximum is on the 2 before last circle, the axis is by gravity center of the max (1 circle around the max)'
        'if the maximum is inside, the axis is by gravity center of the max (2 circles around the max)'
        b=self.initial_values.translation_param
        event_diapason=self.diapason[0]       
        x0_max,y0_max=self.x_mos[event_diapason.sum().idxmax()],self.y_mos[event_diapason.sum().idxmax()] #on mosaic, mm
        d0_from_center=np.sqrt(x0_max**2+y0_max**2)
        if d0_from_center>self.initial_values.d0_center:
            x0, y0 = x0_max, y0_max
            return x0*-b, y0*-b, d0_from_center
        else:
            d_pmt_to_d_max=np.sqrt( (np.array(self.x_mos)-x0_max)**2  + (np.array(self.y_mos)-y0_max)**2  ) #distance of each PMT from the max PMT, m
            if d0_from_center>self.initial_values.d0_before_c and d0_from_center<self.initial_values.d0_center: #if the max PMT is on the circle before the last one
                max_d=self.initial_values.d_six #the minimum distance of the neighors around           
            elif d0_from_center<self.initial_values.d0_before_c: #if the max PMT is inside
                max_d=self.initial_values.d_double_six #the minimum distance of the neighors around        
            number_photons=list(event_diapason[event_diapason>0].sum()) #max of each PMT          
            x_circle_around_max=self.x_mos[np.where(d_pmt_to_d_max<max_d)[0]] #x around the max PMT
            y_circle_around_max=self.y_mos[np.where(d_pmt_to_d_max<max_d )[0]] #y around the max PMT
            number_photons_around_max=np.array(number_photons)[np.where(d_pmt_to_d_max<max_d )[0]] #max PMT around the max PMT           
            x0=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max) #x by gravity center
            y0=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max) #y by gravity center
            return x0*-b, y0*-b, d0_from_center
    
    def length(self): 
        'Return the length of the signal in ns'
        if self.type_analyse=='only_sig': #just the signal without background
            event=self.event #read_mosaic_hits_file()
            event_cells_sup_0=pd.DataFrame({'event':np.array(event.max(axis=1)),'cells':self.cells})
            event_cells_sup_0=event_cells_sup_0[event_cells_sup_0['event']>0]
            cells_sup0=event_cells_sup_0['cells'] #Delete 5% on the right and the left
            lim_inf,lim_sup = np.percentile(np.array(cells_sup0), 5), np.percentile(np.array(cells_sup0), 95) #percent of the impulse saved
            idx_95=self.cells[cells_sup0[cells_sup0>lim_sup].index[0]-1]
            idx_5=self.cells[cells_sup0[cells_sup0<lim_inf].index[-1]+1]
            lenght_impulse=idx_95-idx_5 #lenght
        elif self.type_analyse!='only_sig': #signal with background
            cells_diapason=self.diapason[1]
            lenght_impulse=cells_diapason[-1]-cells_diapason[0]-(self.margin_left+self.margin_right)*self.initial_values.bins  #lenght
        return lenght_impulse   
    
    def amplification(self): 
        'Return the amplified signal in each channel (N, t)'
        find_diapason=self.diapason
        event_diapason,cells_diapason=find_diapason[0],find_diapason[1]
        amplif_signal=event_diapason.rolling(window=self.initial_values.window).sum().dropna()#amplification of the signal by sliding window
        amplif_signal.index=np.arange(0,len(amplif_signal)) #reindex
        amplif_cells=cells_diapason[:-3] #t, ns
        return amplif_signal, amplif_cells #N photons, time  
    
    def noise_pmt(self): #Noise in each PMT
        'Return the noise in piedestal for each channel'
        noise_window=self.event[1:self.pied].rolling(window=self.initial_values.window).sum().dropna() #amplify the piedestal
        if self.type_analyse=='only_sig':
            mean_std=[0]*self.total_number_of_pmt #threshold <bg>+3 std
        else:
            mean_std=noise_window[noise_window>0].mean()+self.initial_values.N_std_background_each*noise_window[noise_window>0].std() #threshold <bg>+3std
        return mean_std

    def front_up_noise(self): #Selection CRITERIA PMT: noise / Front [x,y,t,N] above the noise on snow
        'First criteria for PMT'
        'Return the PMT with N photons above the background'
        N_t_amplification=self.amplification #amplification of the new diapason
        all_window,cells_window=N_t_amplification[0],N_t_amplification[1] 
        noise_pmt=self.noise_pmt
        translation_snow_mos_results=self.front_snow
        x_snow,y_snow,t_path=translation_snow_mos_results[0],translation_snow_mos_results[1],translation_snow_mos_results[2]
        x_y_t_N=pd.DataFrame(
            {'PMT': self.initial_values.num_pmt,'x':x_snow,'y':y_snow,'N_max':list(all_window.max()),
              't_max':cells_window[all_window.idxmax()],'index':all_window.idxmax(),
              'noise_thrs':noise_pmt}).sort_values(by='N_max', ascending=False)
        x_y_t_N=x_y_t_N[x_y_t_N['N_max']>x_y_t_N['noise_thrs']] #save the PMT up to max noise          
        x_y_t_N['t_max']=x_y_t_N['t_max']-t_path[np.array(x_y_t_N['PMT']).astype(int)]        
        return x_y_t_N       
        
    def DFS(self): #Selection CRITERIA PMT: diapason / DFS method: PMT with signal, Front [x,y,t,N] 
        'Second criteria for PMT, Front [x,y,t,N] '
        'Return the PMT that pass through DFS'
        'The impulse time is within the limit of time t ± n of a neighbor '
        x_y_t_N=self.front_up_noise #x,y,t,N        
        x_y_t_N_good=x_y_t_N.iloc[0].to_frame().T #x,y,t,N: axis --> good PMT
        x_y_t_N_maybe=x_y_t_N.iloc[1:] #Others --> don't know        
        PMT_around=self.initial_values.circle_of_pmt_around_pmt.loc[x_y_t_N_good.index[0]].dropna()[1:]  #PMT around the max, delete the central PMT
        #Move them from all to good table
        x_y_t_N_good=pd.concat([x_y_t_N_good,x_y_t_N_maybe.loc[PMT_around[PMT_around.isin(x_y_t_N_maybe['PMT'])]]]) 
        x_y_t_N_maybe=x_y_t_N_maybe.drop(PMT_around[PMT_around.isin(x_y_t_N_maybe['PMT'])]) #Delete 6 neighbors around the max PMT
        for k in range(0,20):
            good_pmt=[]
            bad_pmt=[]
            for i in range(0,len(x_y_t_N_maybe)):
                looked=x_y_t_N_maybe.iloc[i] #number of the PMT in the table 2
                PMT_around=self.initial_values.circle_of_pmt_around_pmt.loc[looked.name].dropna()[1:]  #PMT around the PMT, delete the central PMT
                PMT_around_in_table_1=PMT_around[PMT_around.isin(x_y_t_N_good['PMT'])]                
                if len(PMT_around_in_table_1)==0: #No neighbours of the PMT of table 2 in the table 1
                    continue
                else:
                    mean_time=x_y_t_N_good.loc[PMT_around_in_table_1]['t_max'].mean() #mean of the sure PMT around the examinated PMT                    
                    if looked['t_max'] <= mean_time + self.initial_values.lim_search and looked['t_max'] >= mean_time - self.initial_values.lim_search: 
                        good_pmt.append(looked['PMT'])
                    else:
                        bad_pmt.append(looked['PMT'])
                        continue
            x_y_t_N_good=pd.concat([x_y_t_N_good,x_y_t_N_maybe.loc[good_pmt]])
            x_y_t_N_maybe=x_y_t_N_maybe.drop(good_pmt+bad_pmt) #Delete sure PMT        
        x_y_t_N_good['t_max_on_min']=np.array(x_y_t_N_good['t_max'])-min(np.array(x_y_t_N_good['t_max'])) #Readjust ti on t min
        x_y_t_N_good=x_y_t_N_good.sort_values(by='N_max')       
        return x_y_t_N_good        
        
    def neighbors(self): 
        'Third criteria for PMT, Front [x,y,t,N] '
        'Return the PMT with at least 2 neighbors PMT'
        x_y_t_N_good = self.DFS
        bad_PMT=[]
        for i in range(0,len(x_y_t_N_good)):
            if len(x_y_t_N_good['PMT'][x_y_t_N_good['PMT'].isin(self.initial_values.circle_of_pmt_around_pmt.loc[x_y_t_N_good.iloc[i][0]].dropna())].drop(x_y_t_N_good.iloc[i][0]))<2:
                bad_PMT.append(int(x_y_t_N_good.iloc[i][0]))       
        x_y_t_N_neighbors=x_y_t_N_good.drop(bad_PMT)       
        return x_y_t_N_neighbors
        
    def angles(self):
        'Return theta (rad), phi (rad), a0, a1, a2 '
        'Using multistart to avoid the local minimum'
        'The algorithm is based on the least squares method'
        theta_initt=[0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]*8 #multistart, in rad
        phi_initt=np.arange(0.5,6.5,0.5).tolist()*8       #multistart, in rad
        x_y_t_N_neighbors=self.neighbors     
        x_front,y_front,t_fnc=(np.array(x_y_t_N_neighbors['x']),
                                np.array(x_y_t_N_neighbors['y']),
                                np.array(x_y_t_N_neighbors['t_max_on_min']))    
        N_photons = np.array(x_y_t_N_neighbors['N_max']) ## extract values of the x ,y, t, N of the front on snow           
        def fnc(theta,phi,a0,a1,a2):#определение функции с не заданными параметрами
                x_casc=((np.cos(theta)*np.cos(phi))*(x_front)+(np.cos(theta)*np.sin(phi))*(y_front)) #координаты x фотонов в системе ливня
                y_casc=(-np.sin(phi)*(x_front)+np.cos(phi)*(y_front)) #координаты y фотонов в системе ливня
                z_casc=((np.sin(theta)*np.cos(phi))*(x_front)+(np.sin(theta)*np.sin(phi))*(y_front)) #координаты z фотонов в системе ливня
                R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
                tau=a0+a1*R+a2*R**2+z_casc/self.initial_values.c_ns #аппроксимированое время
                s=(((t_fnc-tau)**2))*N_photons #функцию которую надо аппроксимировать
                Smin=np.sum(s)
                return Smin           
        theta_multi_start=[]
        phi_multi_start=[]       
        a0_multi_start=[]
        a1_multi_start=[]
        a2_multi_start=[]          
        min_s=[]
        for th in range(0,len(theta_initt)):
                param_fix=Minuit(fnc, theta=theta_initt[th],phi=phi_initt[th],a0=0,a1=0,a2=0)
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
        return theta, phi, a0, a1, a2
        
       
        
       
        