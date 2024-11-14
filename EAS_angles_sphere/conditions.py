import pandas as pd
import numpy as np

class Conditions:
    'This class has all conditions (criteria) of events selection'
    def __init__(self,algorithm,initial_values):
        self.algorithm=algorithm
        self.initial_values=initial_values
        self.type_analyse=self.initial_values.type_analyse
        self.intensity_signal=self.algorithm.intensity_signal
        self.algorithm_condition=self.algorithm_condition()
        
    def algorithm_condition(self):
        
        if self.type_analyse!='only_sig' and self.intensity_signal[1]<self.initial_values.condition_intensity:
            a = 'too low signal'
            return
        
        elif self.type_analyse=='only_sig' and self.intensity_signal[0]<self.initial_values.condition_max_impulse:
            a = 'too low signal'
            return
        
        if self.algorithm.diapason[1][0]+self.initial_values.margin_left*self.initial_values.bins<self.initial_values.pied*self.initial_values.bins:
            a = 'event in piedestal'
            return
        
        if self.algorithm.axis[2]>self.initial_values.d0_center:
            a = 'axis too far'
            return
        
        if self.algorithm.lenght_impulse<self.initial_values.condition_length_min or self.algorithm.lenght_impulse>self.initial_values.condition_length_max:
            a = 'lenght too short or long'
            return
        
        elif self.intensity_signal[1]>self.initial_values.condition_calibration or self.algorithm.lenght_impulse>self.initial_values.condition_length_max:
            a = 'calibration event'
            return
        
        if len(self.algorithm.front_up_noise) < self.initial_values.condition_N_pmt:
           a = 'not enough PMT'
           return
         
        if self.algorithm.angles[4]<-0.005:
            a = 'a2<0'
            return
        
        a = 'done'
        
        return a
       