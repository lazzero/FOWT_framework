# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 11:11:01 2021

@author: Giancarlo&Guido
Analyse time series in the frequency domain

## function name: 
    "myPSD"
## inputs: 
    "data" (time series of data to be analysed, e.g. yaw angle in "deg", numpy vector [Nx1])
    "time" (time series of time instants corresponding to the time series of data, in "s", numpy vector [Nx1])
    "padded_length" (number of points to be considered in the fourier transform, integer)
    "folder_path" (path to the folder where the plot files would be saved, if "plot_flag" is set to True)
    "plot_flag" (flag to save plots of the frequency domain analyses)
## outputs:
    "freqs" (frequencies of the frequency domain analysis performed)
    "power_spectrum" (power spectral density)
    "fourier_transform" (fourier transform of the time series         
            If `padded_length` is even, the length of the transformed axis is ``(padded_length/2)+1``.
            If `padded_length` is odd, the length is ``(padded_length+1)/2``.)
    "fourier_transform_norm"
    
"""
import matplotlib.pyplot as plt
import numpy as np

def freq_domain_data(data,time,padded_length=None,folder_path='unused',plot_flag=True,plot_name = 'unused'):
            
    fs = 1/(time[1]-time[0]) # frequency of the sampled data (calculated from the first two instants) [hz]
    
    df = 1/max(time) # minimum resolution in frequency of the sampled data (calculated from last instant) [hz]
    
    N = len(data) # number of sampled data
    
    if padded_length==None:
        padded_length = N
    elif isinstance(padded_length,int):
        pass
    else: 
        raise("Non integer number of points") 
    
    data = data-np.mean(data)
    
    fourier_transform = np.fft.rfft(data,n=padded_length) # one-dimensional discrete Fourier Transform for real input (sampled data)
    
    # define length of the trasformed data "N_trasf" based on odd/even padded length
    
    if padded_length % 2 == 0:
        N_trasf = int((padded_length/2)+1)
    else:
        N_trasf = int((padded_length+1)/2)
    
    # calculate power spectrum and frequency vector
    
    fourier_transform_squared=(np.abs(fourier_transform))**2
    
    power_spectrum = (fourier_transform_squared/(N)**2)/df

    freqs = np.linspace(0, fs/2, N_trasf) # vector of the frequencies analysed
    
    fourier_transform_norm = np.abs(fourier_transform)/(N/2) # https://dsp.stackexchange.com/questions/66058/am-i-supposed-to-normalize-fft-in-python
    
    freq_PSDpeak = freqs[np.argmax(power_spectrum)]
    PSDpeak = max(power_spectrum)

    freq_FFTpeak = freqs[np.argmax(fourier_transform_norm)]
    FFTpeak = max(fourier_transform_norm)    
    # plot frequency domain analyses data
    
    if plot_flag:
        plt.xlabel('f (Hz)')
        plt.ylabel('PSD (sig^2/Hz)')
        plt.xlim([0,0.25])
        plt.title('{name} PSD'.format(name=plot_name))
        plt.plot(freqs,power_spectrum)
        plt.text(freq_PSDpeak+freq_PSDpeak*0.125,PSDpeak-PSDpeak*0.05,f"f={freq_PSDpeak:5.3f} Hz")
        plt.text(freq_PSDpeak+freq_PSDpeak*0.125,PSDpeak-PSDpeak*0.10,f"T={1/freq_PSDpeak:5.2f} s")
        plt.savefig(folder_path+"\\output_{name}_PSD.png".format(name=plot_name)) #save as png
        plt.clf()
        
        plt.xlabel('f (Hz)')
        plt.ylabel('FFT (sig)')
        plt.xlim([0,0.25])
        plt.title('{name} FFT Normalized'.format(name=plot_name))
        plt.plot(freqs,fourier_transform_norm)
        plt.text(freq_FFTpeak+freq_FFTpeak*0.125,FFTpeak-FFTpeak*0.05,f"f={freq_FFTpeak:5.3f} Hz")
        plt.text(freq_FFTpeak+freq_FFTpeak*0.125,FFTpeak-FFTpeak*0.10,f"T={1/freq_FFTpeak:5.2f} s")
        plt.savefig(folder_path+"\\output_{name}_fft.png".format(name=plot_name)) #save as png
        plt.clf()
        
    
    return freqs, power_spectrum , fourier_transform, fourier_transform_norm