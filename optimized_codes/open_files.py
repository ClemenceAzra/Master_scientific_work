import numpy as np
import pandas as pd


class Openfiles:
    
    'This class opens different types of files and output DataFrame massive with the size: response*N PMT'
    'each value corresponds to the number of photons '
    'made for the SPHERE-2/ SPHERE-3 telescopes and the type of analyse "only signal, electronic, experiment"'
    'you don t have to touch anything expect you want to add another type of file. Then create new "if"'
    
    def __init__(self, initial_values):
        self.initial_values=initial_values
        self.path=self.initial_values.path
        self.type_analyse=self.initial_values.type_analyse
        self.telescope=self.initial_values.telescope
        self.event=self.read_mosaic_hits_file()
    
    def read_mosaic_hits_file(self): 
        q = np.zeros((self.initial_values.total_number_of_pmt, self.initial_values.response))  # empty array 109 PMT * 1020 cells
       
        def create_dataframe(self): #for only signal
            for i in np.arange(self.initial_values.num_pmt.shape[0]):
                impulse = np.array(pd_event)[:,1][np.where(np.array(pd_event)[:,0]==i)]
                n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
                n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
                q[i]=n_phot_shift
            return q
        def substract_piedestal(self): #for electronic and experiment
             for pmt in range(event.shape[1]):
                 q.T[0::2, pmt] = event[0::2, pmt] - np.mean(event[0:self.initial_values.pied:2, pmt])
                 q.T[1::2, pmt] = event[1::2, pmt] - np.mean(event[1:self.initial_values.pied:2, pmt])
             for i in range(event.shape[1]):
                 q.T[i]=q.T[i]*self.initial_values.coeff_amp[1].iloc[i]
             return q
    
        if self.initial_values.telescope =='SPHERE-2':      
            if self.type_analyse=='only_sig':
                shift = 500 #first bin of the signal
                cells_here=np.arange(0, 1021*self.initial_values.bins, self.initial_values.bins)
                a = pd.read_csv(self.path, header=None, skiprows=[0],sep='\s+')
                pd_event = np.array(pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[self.initial_values.num_t_event]-min(a[self.initial_values.num_t_event]))})) #creation of DataFrame to facilite calculations
                return pd.DataFrame(create_dataframe(self).T)

            elif self.type_analyse=='electronic':
                event = pd.read_csv(self.path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip').to_numpy()
                return pd.DataFrame(substract_piedestal(self).T)

            elif self.type_analyse=='experiment':
                event = pd.read_csv(self.path, skiprows=40, skipfooter=2144,engine='python',sep='\s+',header=None)
                event = event.drop(columns=[0,110,111,112])
                event.columns -= 1 #name of PMT from 0 to 108
                event = event.to_numpy()
                return pd.DataFrame(substract_piedestal(self).T)

            
            
            
            
            
            
                