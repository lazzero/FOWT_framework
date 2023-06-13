# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 12:46:19 2022

@author: guila
"""
import os
import pandas as pd
from timetofreqdomain import freq_domain_data
from scipy.signal import find_peaks

file = "C:/FOWT_optim_test/sims/template_input_files/TimeMomentSeries.dat"

data = pd.read_csv(file,'\s',skiprows=6)

time_col = data['#'][1:].to_numpy(na_value=0)
data_col = data['(N-m).1'][1:].to_numpy(na_value=0)

[f,psd,fft,fft_n] = freq_domain_data(data_col,time_col,folder_path=os.getcwd(),plot_flag = True)

fft_n_max = max(fft_n)
min_peak_height = fft_n_max-fft_n_max*0.01

[peak_index,properties] = find_peaks(fft_n,height=(min_peak_height))

f_exc = f[peak_index] # excitation frequency in Hz
T_exc = 1/f_exc # excitation frequency in s

print('\nexcitation frequency: %.4f Hz \n' % (f_exc))
print('excitation period: %.1f s \n'% (T_exc))