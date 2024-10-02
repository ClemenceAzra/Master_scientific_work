import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class Results:
    'This class return the results: error, figures, and save them'
    
    def __init__(self, initial_values, algorithm):
        self.algorithm=algorithm
        self.initial_values=initial_values
        self.angles=self.algorithm.angles
        self.delta=self.delta()
        self.fig_sum_impulse=self.fig_sum_impulse()
        self.front=self.front()
        
        
    def delta(self): 
       'Return the error of angle determination - delta'
       if self.initial_values.type_analyse!='experiment':
           theta_fit=self.angles[0]
           phi_fit=self.angles[1]
           theta_real=self.initial_values.real_theta
           phi_real=self.initial_values.real_phi
           cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
           cy = np.sin(theta_real) * np.sin(phi_real)
           cz = np.cos(theta_real)
           cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
           cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
           cz_fit = np.cos(theta_fit)
           delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
           return delta
       if self.initial_values.type_analyse=='experiment':
           return np.nan


    def save_results(self):
        'Save results: angles, axis, length...'
        angles_and_param=self.angles
        axis_values=self.axis
        lenght=self.algorithm.lenght
        theta, phi, a0, a1, a2 = angles_and_param[0], angles_and_param[1], angles_and_param[2], angles_and_param[3], angles_and_param[4]
        axis_x, axis_y = axis_values[0], axis_values[1]
        table_final=pd.DataFrame({'theta':theta, 'phi':phi,'a0':a0,'a1':a1,'a2':a2, 'axis_x':axis_x, 'axis_y':axis_y,'length':lenght})
        # table_final.to_csv('/Users/clemence/Documents/results)
    
    def fig_sum_impulse(self):
        'Draw the summed impulse in each PMT'
        event = self.algorithm.event
        noise_mean_std=self.algorithm.diapason[2]
        
        plt.bar(self.initial_values.cells,event.sum(axis=1),width=50)
        if self.initial_values.type_analyse!='only_sig':
            plt.axhline(y=noise_mean_std,color='darkorange',linewidth=3)
        plt.grid()
        plt.ylabel('∑ Ni фот',fontsize=14)
        plt.xlabel('t, нс',fontsize=14)
        plt.xlim(0,self.initial_values.cells[-1])
        plt.ylim(min(event.sum(axis=1))-20,max(event.sum(axis=1))+20)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.show() 
        
    def front(self):      
        'Draw the 3D front on snow'
        front_3d = self.algorithm.neighbors
        fig = plt.figure(figsize=(8,10))
        ax = plt.axes(projection="3d")
        x_1=front_3d['x']
        y_1=front_3d['y']
        t_1=front_3d['t_max_on_min']
        ax.scatter(x_1,y_1,t_1,s=10,label='Experimental front',color='blue')
        # # ax.scatter(x_1,y_1,tau_fit,s=10,label='Аппроксимация',color='red')
        ax.view_init(10,140)  #For diploma: sph2 (10.80), sph3 (10,115)
        ax.set_xlabel("x, м",fontsize=14)
        ax.set_ylabel("y, м",fontsize=14)
        ax.set_zlabel("t, нс",fontsize=14)
        ax.xaxis.labelpad = 10
        ax.yaxis.labelpad = 10
        ax.zaxis.labelpad = 5
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(np.arange(-200,300,100))
        plt.yticks(np.arange(-200,300,100))
        ax.set_box_aspect(aspect=None, zoom=0.8)
        # # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/3d_front_sph2.pdf',bbox_inches='tight')
        plt.show() #3d_front_sph2
        
        
        
        
        
        