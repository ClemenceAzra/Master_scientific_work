import pandas as pd
import numpy as np

class Initial_values:
    'This class inputs initial values linked to the telescope geometry, files... '
    'You SHOULD change the path of files, parameter <path_init_values>' 
    'You can change others constant parameters if necessary'
    
    def __init__(self, type_analyse, telescope, dirname_full,nuclei,En,H):
        self.type_analyse = type_analyse
        self.telescope = telescope
        self.altitude=H
        self.pied=400
        self.bins = 12.5     #cell, ns
        self.response = 1020
        self.c_ns=0.3
        self.path = dirname_full
        
        path_init_values='/Users/clemence/Documents/Магистратура_наука/Научная_работа/Data/4_data_init' #YOU CAN CHANGE the path of the initial values (PMT coordinates, etc...)

        
        if self.telescope =='SPHERE-2':         ### for ! SPHERE-2 !
            self.total_number_of_pmt = 109              ##  total number of PMT
            self.num_t_event = 4                        ##  column of time in initial file (only signal)
            self.d0_center=192                   #last circle on mosaic
            self.d0_before_c=154                 #before last circle on mosaic
            self.d_six=50                        #distance PMT and 1st circle around
            self.d_double_six=100                #distance PMT and 2nd circles around
            self.lim_search=100                  #diapason for search the impulse by DFS, ns
            self.translation_param = 0.9051146902370498  #parameter for translation mossaic --> snow

            self.coeff_amp = pd.read_csv(f'{path_init_values}/coeff_amp.csv', header=None, sep='\s+')          #Coeff. amplification
            x_y_mos=pd.read_csv(f'{path_init_values}/mosaic_pmt_coords_df.csv')  
            self.x_mos=np.array(x_y_mos['x'])    #x on mosaic, mm
            self.y_mos=np.array(x_y_mos['y'])    #y on mosaic, mm
            self.circle_of_pmt_around_pmt=pd.read_csv(f'{path_init_values}/circle_of_pmt_around_pmt.csv')  #Neighbors of PMT
            self.real_theta=np.array(pd.read_csv(f'{path_init_values}/angles_theta',header=None,sep='\s+')[0])[0]
            self.real_phi=np.array(pd.read_csv(f'{path_init_values}/angles_phi',header=None,sep='\s+')[0])[0]

            if self.type_analyse == 'experiment':
                self.pied = 395
                self.response = 915
            
        if self.telescope=='SPHERE-3':           ### for ! SPHERE-3 !
            self.total_number_of_pmt = 379              ##  total number of PMT
            self.num_t_event = 5                        ##  column of time in initial file (only signal)
            self.d0_center=292                   #last circle on mosaic
            self.d0_before_c=262                 #before last circle on mosaic
            self.d_six=32                        #distance PMT and 1st circle around
            self.d_double_six=70                 #distance PMT and 2nd circles around
            self.lim_search=50                   #diapason for search the impulse by DFS, ns
            self.translation_param = 0.5657087664485793  #parameter for translation mossaic --> snow

            pixel_data = pd.read_csv(f'{path_init_values}/SPHERE3_pixel_data_A.dat', header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+')  #Coord. PMT on mosaic, mm
            pixel_data['segment'] = pixel_data.index // 7 #== add num segment
            self.x_mos=np.array(pixel_data.groupby('segment').mean()['x']) #x on mosaic
            self.y_mos=np.array(pixel_data.groupby('segment').mean()['y']) #y on mosaic
            self.circle_of_pmt_around_pmt=pd.read_csv(f'{path_init_values}/circle_of_pmt_around_pmt_sph3.csv')  #Neighbors of PMT  
            self.real_phi=np.array(pd.read_csv(f'{path_init_values}/angels_phi_{self.nuclei}_sph3',header=None,sep='\s+')[0])
            self.real_theta=15    

            if self.type_analyse=='electronic':  
                self.pied=100                    #piedestal
                self.response=500                #lenght of the response
                
        self.num_pmt = np.arange(0, self.total_number_of_pmt) #number of PMT for 0 to the last        
        self.cells = np.arange(0, self.response * self.bins, self.bins)  #length of the electronic response, cells*response
         
        self.condition_intensity=2
        self.condition_max_impulse=50
        self.condition_length_max=2000
        self.condition_length_min=100
        self.condition_calibration=10000
        self.condition_N_pmt=1

            