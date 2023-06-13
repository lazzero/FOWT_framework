# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 11:11:01 2021

@author: Giancarlo
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.signal import find_peaks

def myPSD(data,time,padded_length=None):
    
    sampling_rate = 1/(time[1]-time[0])
    
    df=1/max(time)
    
    if padded_length==None:
        padded_length=len(data)
    elif isinstance(padded_length,int):
        pass 
    else: 
        raise("Non integer number of points")
    
    
    fourier_transform = np.fft.rfft(data,n=padded_length)

    #fourier_transform_squared = fourier_transform * np.conj(fourier_transform)
    
    fourier_transform_squared=(np.abs(fourier_transform))**2
    
    power_spectrum = (fourier_transform_squared/(len(data))**2)/df

    freqs = np.linspace(0, sampling_rate/2, len(power_spectrum))
    
    return freqs, power_spectrum, fourier_transform

if __name__=="__main__":
    t=np.linspace(0.0,100.0,10000)
    xx=10.0*np.sin(2*np.pi*10*t)+5.0*np.sin(2*np.pi*25*t)
    freqs, PSD, fft=myPSD(xx,t,padded_length=10000)
    
    f2 = plt.figure()
    plt.plot(t,xx)
    f1 = plt.figure()
    plt.plot(freqs,PSD)
    plt.xlabel("f (Hz)")
    plt.ylabel("PSD")
    f_max=freqs[np.argmax(PSD)]
    PSD_max=np.max(PSD)
    plt.text(f_max,PSD_max,"f_max={:.3f} Hz --- PSD_max={:.1f} ".format(f_max,PSD_max))
    
    f3 = plt.figure()
    plt.plot(freqs,np.abs(fft)/len(t)*2)
    plt.xlabel("f (Hz)")
    plt.ylabel("fft")
    f_max1=freqs[np.argmax(np.abs(fft))]
    fft_max=np.max(np.abs(fft)/len(t)*2)
    plt.text(f_max1,fft_max,"f_max={:.3f} Hz --- fft_max={:.1f} ".format(f_max1,fft_max))
    
    df=pd.read_csv('Myaw_Helix_4degrees_8times.dat',sep=' ')
    
    peaks, _ = find_peaks(df.Myaw,width=40)
    peaks_index = np.argsort(peaks)
    plt.plot(df.t[peaks], df.Myaw[peaks], "x")
    t_yaw = np.mean(np.diff(df.t[peaks]))
    plt.plot(df.t,df.Myaw)
