# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 12:20:43 2023

@author: guila
"""

from pyFAST.input_output.fast_output_file import FASTOutputFile
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


#Options
params = {'text.usetex' : True,
          'font.size' : 15,
          'font.family' : 'palatino',
          'figure.figsize' : (10,5)
          }

plt.rcParams.update(params) 

def freq_domain_data(data,time,padded_length=None,folder_path='unused',plot_flag=True, plot_name = 'unused'):
            
    fs = 1/(time[1]-time[0]) # frequency of the sampled data (calculated from the first two instants) [hz]
    
    df = 1/max(time) # minimum resolution in frequency of the sampled data (calculated from last instant) [hz]
    
    N = len(data) # number of sampled data
    
    if padded_length==None:
        padded_length = N
    elif isinstance(padded_length,int):
        pass
    else: 
        raise("Non integer number of points") 
    
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
    
    freq_FFTpeak = freqs[np.argmax(fourier_transform_norm)]
    freq_PSDpeak = freqs[np.argmax(power_spectrum)]

    FFTpeak = max(fourier_transform_norm)
    PSDpeak = max(power_spectrum)
    # plot frequency domain analyses data
    
    if plot_flag:
        plt.xlabel('$f (Hz)$')
        plt.ylabel('$PSD (sig^2/Hz)$')
        plt.xlim([0,0.25])
        plt.title('{name} PSD'.format(name=plot_name))
        plt.plot(freqs,power_spectrum)
        plt.text(freq_PSDpeak+freq_PSDpeak*0.125,PSDpeak-PSDpeak*0.05,f"f={freq_PSDpeak:5.3f} Hz")
        plt.text(freq_PSDpeak+freq_PSDpeak*0.125,PSDpeak-PSDpeak*0.10,f"T={1/freq_PSDpeak:5.2f} s")
        plt.text(freq_PSDpeak+freq_PSDpeak*0.125,PSDpeak-PSDpeak*0.15,f"max={PSDpeak:5.2f} $(sig^2/Hz)$")
        plt.savefig(folder_path+"\\output_{name}_PSD.png".format(name=plot_name)) #save as png
        plt.clf()
        
        plt.xlabel('f (Hz)')
        plt.ylabel('$FFT (sig)$')
        plt.xlim([0,0.25])
        plt.title('{name} FFT Normalized'.format(name=plot_name))
        plt.plot(freqs,fourier_transform_norm)
        plt.text(freq_FFTpeak+freq_FFTpeak*0.125,FFTpeak-FFTpeak*0.05,f"f={freq_FFTpeak:5.3f} Hz")
        plt.text(freq_FFTpeak+freq_FFTpeak*0.125,FFTpeak-FFTpeak*0.10,f"T={1/freq_FFTpeak:5.2f} s")
        plt.text(freq_FFTpeak+freq_FFTpeak*0.125,FFTpeak-FFTpeak*0.15,f"max={FFTpeak:5.2f} (sig)")
        plt.savefig(folder_path+"\\output_{name}_fft.png".format(name=plot_name)) #save as png
        plt.clf()
        
    
    return freqs, power_spectrum , fourier_transform, fourier_transform_norm

folder_path = r"C:\FOWT_optim_SOFTWIND_rev5\FOWT_optim_test\sims\baseline_input_files_triple_spar_heave_thrust"
output_filename = r'\DTU10MW3Spar_Param_MPHC_par.outb'

outdata=FASTOutputFile(folder_path+output_filename).toDataFrame()

evalTime = 0

outdata = outdata[outdata["Time_[s]"]>evalTime]


rdx_mean = np.mean(outdata["PtfmRoll_[deg]"].to_numpy()[:])
rdy_mean = np.mean(outdata["PtfmPitch_[deg]"].to_numpy()[:])
rdz_mean = np.mean(outdata["PtfmYaw_[deg]"].to_numpy()[:])
tdx_mean = np.mean(outdata["PtfmSurge_[m]"].to_numpy()[:])
tdy_mean = np.mean(outdata["PtfmSway_[m]"].to_numpy()[:])
tdz_mean = np.mean(outdata["PtfmHeave_[m]"].to_numpy()[:])

rdx_std = np.std(outdata["PtfmRoll_[deg]"].to_numpy()[:])
rdy_std = np.std(outdata["PtfmPitch_[deg]"].to_numpy()[:])
rdz_std = np.std(outdata["PtfmYaw_[deg]"].to_numpy()[:])
tdx_std = np.std(outdata["PtfmSurge_[m]"].to_numpy()[:])
tdy_std = np.std(outdata["PtfmSway_[m]"].to_numpy()[:])
tdz_std = np.std(outdata["PtfmHeave_[m]"].to_numpy()[:])

rdx_max = np.abs(max(outdata["PtfmRoll_[deg]"].to_numpy()[:]))
rdy_max = np.abs(max(outdata["PtfmPitch_[deg]"].to_numpy()[:]))
rdz_max = np.abs(max(outdata["PtfmYaw_[deg]"].to_numpy()[:]))
tdx_max = np.abs(max(outdata["PtfmSurge_[m]"].to_numpy()[:]))
tdy_max = np.abs(max(outdata["PtfmSway_[m]"].to_numpy()[:]))
tdz_max = np.abs(max(outdata["PtfmHeave_[m]"].to_numpy()[:]))

time=outdata["Time_[s]"].to_numpy()[:]

roll=outdata["PtfmRoll_[deg]"].to_numpy()[:]
pitch=outdata["PtfmPitch_[deg]"].to_numpy()[:]
yaw=outdata["PtfmYaw_[deg]"].to_numpy()[:]
surge=outdata["PtfmSurge_[m]"].to_numpy()[:]
sway=outdata["PtfmSway_[m]"].to_numpy()[:]
heave=outdata["PtfmHeave_[m]"].to_numpy()[:]

rdx=outdata["PtfmRoll_[deg]"].to_numpy()[:]-rdx_mean
rdy=outdata["PtfmPitch_[deg]"].to_numpy()[:]-rdy_mean
rdz=outdata["PtfmYaw_[deg]"].to_numpy()[:]-rdz_mean
tdx=outdata["PtfmSurge_[m]"].to_numpy()[:]-tdx_mean
tdy=outdata["PtfmSway_[m]"].to_numpy()[:]-tdy_mean
tdz=outdata["PtfmHeave_[m]"].to_numpy()[:]-tdz_mean

thrust = outdata["RotThrust_[kN]"].to_numpy()[:]
#power = outdata["RotPwr_[kW]"].to_numpy()[:]
rotspeed = outdata["RotSpeed_[rpm]"].to_numpy()[:]

thrust_mean = np.mean(outdata["RotThrust_[kN]"].to_numpy()[:])
#power_mean = np.mean(outdata["RotPwr_[kW]"].to_numpy()[:])
rotspeed_mean = np.mean(outdata["RotSpeed_[rpm]"].to_numpy()[:])

thrust_std = np.std(outdata["RotThrust_[kN]"].to_numpy()[:])
#power_std = np.std(outdata["RotPwr_[kW]"].to_numpy()[:])
rotspeed_std = np.std(outdata["RotSpeed_[rpm]"].to_numpy()[:])

thrust_max = np.abs(max(outdata["RotThrust_[kN]"].to_numpy()[:]))
#power_max = np.abs(max(outdata["RotPwr_[kW]"].to_numpy()[:]))
rotspeed_max = np.abs(max(outdata["RotSpeed_[rpm]"].to_numpy()[:]))

# Frequency domain analysis
freqs,Pxx1,fft1,fft1_normalized = freq_domain_data(rdx,time, None,folder_path=folder_path,plot_flag=True,plot_name ='roll')

# Function evaluation
f_max=freqs[np.argmax(fft1_normalized)]

# Frequency domain analysis
freqs,Pxx1,fft1,fft1_normalized = freq_domain_data(rdy,time, None,folder_path=folder_path,plot_flag=True,plot_name ='pitch')

# Function evaluation
f_max=freqs[np.argmax(fft1_normalized)]

# Frequency domain analysis
freqs,Pxx1,fft1,fft1_normalized = freq_domain_data(rdz,time, None,folder_path=folder_path,plot_flag=True,plot_name ='yaw')

# Function evaluation
f_max=freqs[np.argmax(fft1_normalized)]

# Frequency domain analysis
freqs,Pxx1,fft1,fft1_normalized = freq_domain_data(tdx,time, None,folder_path=folder_path,plot_flag=True,plot_name ='surge')

# Function evaluation
f_max=freqs[np.argmax(fft1_normalized)]

# Frequency domain analysis
freqs,Pxx1,fft1,fft1_normalized = freq_domain_data(tdy,time, None,folder_path=folder_path,plot_flag=True,plot_name ='sway')

# Function evaluation
f_max=freqs[np.argmax(fft1_normalized)]

# Frequency domain analysis
freqs,Pxx1,fft1,fft1_normalized = freq_domain_data(tdz,time, None,folder_path=folder_path,plot_flag=True,plot_name ='heave')

# Function evaluation
f_max=freqs[np.argmax(fft1_normalized)]


# Plot time histories 

plt.xlabel('time (s)')
plt.ylabel('roll (deg)')
plt.title('Roll Time History')
plt.plot(time,roll)
t = plt.text((time[0]+time[-1])/2,rdx_mean,f"mean = {rdx_mean:5.3f} deg")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(rdx_mean+rdx_max)/2,f"std = {rdx_std:5.3f} deg")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,rdx_max,f"max = {rdx_max:5.3f} deg")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_roll_time_history.png") #save as png
plt.clf()

plt.xlabel('time (s)')
plt.ylabel('pitch (deg)')
plt.title('Pitch Time History')
plt.plot(time,pitch)
t = plt.text((time[0]+time[-1])/2,rdy_mean,f"mean = {rdy_mean:5.3f} deg")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(rdy_mean+rdy_max)/2,f"std = {rdy_std:5.3f} deg")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,rdy_max,f"max = {rdy_max:5.3f} deg")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_pitch_time_history.png") #save as png
plt.clf()

plt.xlabel('time (s)')
plt.ylabel('yaw (deg)')
plt.title('Yaw Time History')
plt.plot(time,yaw)
t = plt.text((time[0]+time[-1])/2,rdz_mean,f"mean = {rdz_mean:5.3f} deg")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(rdz_mean+rdz_max)/2,f"std = {rdz_std:5.3f} deg")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,rdz_max,f"max = {rdz_max:5.3f} deg")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_yaw_time_history.png") #save as png
plt.clf()


plt.xlabel('time (s)')
plt.ylabel('surge (m)')
plt.title('Surge Time History')
plt.plot(time,surge)
t = plt.text((time[0]+time[-1])/2,tdx_mean,f"mean = {tdx_mean:5.3f} m")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(tdx_mean+tdx_max)/2,f"std = {tdx_std:5.3f} m")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,tdx_max,f"max = {tdx_max:5.3f} m")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_surge_time_history.png") #save as png
plt.clf()

plt.xlabel('time (s)')
plt.ylabel('sway (m)')
plt.title('Sway Time History')
plt.plot(time,sway)
t = plt.text((time[0]+time[-1])/2,tdy_mean,f"mean = {tdy_mean:5.3f} m")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(tdy_mean+tdy_max)/2,f"std = {tdy_std:5.3f} m")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,tdy_max,f"max = {tdy_max:5.3f} m")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_sway_time_history.png") #save as png
plt.clf()

plt.xlabel('time (s)')
plt.ylabel('heave (m)')
plt.title('Heave Time History')
plt.plot(time,heave)
t = plt.text((time[0]+time[-1])/2,tdz_mean,f"mean = {tdz_mean:5.3f} m")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(tdz_mean+tdz_max)/2,f"std = {tdz_std:5.3f} m")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,tdz_max,f"max = {tdz_max:5.3f} m")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_heave_time_history.png") #save as png
plt.clf()

plt.xlabel('time (s)')
plt.ylabel('thrust (kN)')
plt.title('Thrust Time History')
plt.plot(time,thrust)
t = plt.text((time[0]+time[-1])/2,thrust_mean,f"mean = {thrust_mean:5.3f} kN")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(thrust_mean+thrust_max)/2,f"std = {thrust_std:5.3f} kN")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,thrust_max,f"max = {thrust_max:5.3f} kN")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_thrust_time_history.png") #save as png
plt.clf()

#plt.xlabel('time (s)')
#plt.ylabel('power (kW)')
#plt.title('Power Time History')
#plt.plot(time,power)
#t = plt.text((time[0]+time[-1])/2,power_mean,f"mean = {power_mean:5.3f} kW")
#t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
#s = plt.text((time[0]+time[-1])/2,(power_mean+power_max)/2,f"std = {power_std:5.3f} kW")
#s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
#c = plt.text((time[0]+time[-1])/2,power_max,f"max = {power_max:5.3f} kW")
#c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
#plt.xlim(evalTime, time[-1])
#plt.savefig(folder_path+"\\output_power_time_history.png") #save as png
#plt.clf()

plt.xlabel('time (s)')
plt.ylabel('rot. speed (RPM)')
plt.title('Rot. Speed Time History')
plt.plot(time,rotspeed)
t = plt.text((time[0]+time[-1])/2,rotspeed_mean,f"mean = {rotspeed_mean:5.3f} RPM")
t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
s = plt.text((time[0]+time[-1])/2,(rotspeed_mean+rotspeed_max)/2,f"std = {rotspeed_std:5.3f} RPM")
s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
c = plt.text((time[0]+time[-1])/2,rotspeed_max,f"max = {rotspeed_max:5.3f} RPM")
c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
plt.xlim(evalTime, time[-1])
plt.savefig(folder_path+"\\output_rotspeed_time_history.png") #save as png
plt.clf()

bladepitch_output = False

if bladepitch_output == True:

    bladepitch1 = outdata["PtchPMzc1_[deg]"].to_numpy()[:]
    bladepitch2 = outdata["PtchPMzc2_[deg]"].to_numpy()[:]
    bladepitch3 = outdata["PtchPMzc3_[deg]"].to_numpy()[:]

    bladepitch1_mean = np.mean(outdata["PtchPMzc1_[deg]"].to_numpy()[:])
    bladepitch2_mean = np.mean(outdata["PtchPMzc2_[deg]"].to_numpy()[:])
    bladepitch3_mean = np.mean(outdata["PtchPMzc3_[deg]"].to_numpy()[:])

    bladepitch1_std = np.std(outdata["PtchPMzc1_[deg]"].to_numpy()[:])
    bladepitch2_std = np.std(outdata["PtchPMzc2_[deg]"].to_numpy()[:])
    bladepitch3_std = np.std(outdata["PtchPMzc3_[deg]"].to_numpy()[:])

    bladepitch1_max =  np.abs(max(outdata["PtchPMzc1_[deg]"].to_numpy()[:]))
    bladepitch2_max =  np.abs(max(outdata["PtchPMzc2_[deg]"].to_numpy()[:]))
    bladepitch3_max =  np.abs(max(outdata["PtchPMzc3_[deg]"].to_numpy()[:]))

    plt.xlabel('time (s)')
    plt.ylabel('blade 1 pitch (deg)')
    plt.title('Blade 1 pitch')
    plt.plot(time,bladepitch1)
    t = plt.text((time[0]+time[-1])/2,bladepitch1_mean,f"mean = {bladepitch1_mean:5.3f} deg")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    s = plt.text((time[0]+time[-1])/2,(bladepitch1_mean+bladepitch1_max)/2,f"std = {bladepitch1_std:5.3f} deg")
    s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    c = plt.text((time[0]+time[-1])/2,bladepitch1_max,f"max = {bladepitch1_max:5.3f} deg")
    c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_bladepitch1_time_history.png") #save as png
    plt.clf()

    plt.xlabel('time (s)')
    plt.ylabel('blade 2 pitch (deg)')
    plt.title('Blade 2 pitch')
    plt.plot(time,bladepitch2)
    t = plt.text((time[0]+time[-1])/2,bladepitch2_mean,f"mean = {bladepitch2_mean:5.3f} deg")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    s = plt.text((time[0]+time[-1])/2,(bladepitch2_mean+bladepitch2_max)/2,f"std = {bladepitch2_std:5.3f} deg")
    s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    c = plt.text((time[0]+time[-1])/2,bladepitch2_max,f"max = {bladepitch2_max:5.3f} deg")
    c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_bladepitch2_time_history.png") #save as png
    plt.clf()

    plt.xlabel('time (s)')
    plt.ylabel('blade 2 pitch (deg)')
    plt.title('Blade 2 pitch')
    plt.plot(time,bladepitch3)
    t = plt.text((time[0]+time[-1])/2,bladepitch3_mean,f"mean = {bladepitch3_mean:5.3f} deg")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    s = plt.text((time[0]+time[-1])/2,(bladepitch3_mean+bladepitch3_max)/2,f"std = {bladepitch3_std:5.3f} deg")
    s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    c = plt.text((time[0]+time[-1])/2,bladepitch3_max,f"max = {bladepitch3_max:5.3f} deg")
    c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_bladepitch3_time_history.png") #save as png
    plt.clf()

moor_output = False

if moor_output == True:

    folder_path2 = r"C:\Users\Utente\Desktop\Work_Guido\Projects\FLOATECH\WP4\softwind_optimization_helix\optimization_2\OpenFAST_optimum_2_ws10_mom100"
    output_filename2 = r'\SW_SPAR_STATIC_EQUILIBRIUM_EQ3MOORING.MD.out'
    outdata2=FASTOutputFile(folder_path2+output_filename2).toDataFrame()
    outdata2 = outdata2[outdata2["Time_[s]"]>evalTime]
    timeMD=outdata2["Time_[s]"].to_numpy()[:]

    #print(outdata2)

    CON1FZ_mean = np.mean(outdata2["CON1FZ_[N]"])
    CON2FZ_mean = np.mean(outdata2["CON2FZ_[N]"])
    CON3FZ_mean = np.mean(outdata2["CON3FZ_[N]"])
    CON1FZ = outdata2["CON1FZ_[N]"]
    CON2FZ = outdata2["CON2FZ_[N]"]
    CON3FZ = outdata2["CON3FZ_[N]"]

    FAIRTEN1_mean = np.mean(outdata2["FAIRTEN1_[N]"].to_numpy()[:])/1000
    FAIRTEN2_mean = np.mean(outdata2["FAIRTEN2_[N]"].to_numpy()[:])/1000
    FAIRTEN3_mean = np.mean(outdata2["FAIRTEN3_[N]"].to_numpy()[:])/1000
    FAIRTEN1 = outdata2["FAIRTEN1_[N]"].to_numpy()[:]/1000
    FAIRTEN2 = outdata2["FAIRTEN2_[N]"].to_numpy()[:]/1000
    FAIRTEN3 = outdata2["FAIRTEN3_[N]"].to_numpy()[:]/1000
    FAIRTEN1_std = np.std(outdata2["FAIRTEN1_[N]"].to_numpy()[:])/1000
    FAIRTEN2_std = np.std(outdata2["FAIRTEN2_[N]"].to_numpy()[:])/1000
    FAIRTEN3_std = np.std(outdata2["FAIRTEN3_[N]"].to_numpy()[:])/1000
    FAIRTEN1_max = np.abs(max(outdata2["FAIRTEN1_[N]"].to_numpy()[:]))/1000
    FAIRTEN2_max = np.abs(max(outdata2["FAIRTEN2_[N]"].to_numpy()[:]))/1000
    FAIRTEN3_max = np.abs(max(outdata2["FAIRTEN3_[N]"].to_numpy()[:]))/1000


    # Plot time histories 
        # 
    plt.xlabel('time (s)')
    plt.ylabel('CON1FZ (N)')
    plt.title('CON1FZ Time History')
    plt.plot(timeMD,CON1FZ)
    plt.plot([timeMD[0],timeMD[-1]],[CON1FZ_mean,CON1FZ_mean],color="r")
    t = plt.text((timeMD[0]+timeMD[-1])/2,CON1FZ_mean,f"mean = {CON1FZ_mean:5.3f} N")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, timeMD[-1])
    plt.savefig(folder_path+"\\output_CON1FZ_time_history.png") #save as png
    plt.clf()

    plt.xlabel('time (s)')
    plt.ylabel('CON2FZ (N)')
    plt.title('CON2FZ Time History')
    plt.plot(timeMD,CON2FZ)
    plt.plot([timeMD[0],timeMD[-1]],[CON2FZ_mean,CON2FZ_mean],color="r")
    t = plt.text((timeMD[0]+timeMD[-1])/2,CON2FZ_mean,f"mean = {CON2FZ_mean:5.3f} N")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, timeMD[-1])
    plt.savefig(folder_path+"\\output_CON2FZ_time_history.png") #save as png
    plt.clf()

    plt.xlabel('time (s)')
    plt.ylabel('CON3FZ (N)')
    plt.title('CON3FZ Time History')
    plt.plot(timeMD,CON3FZ)
    plt.plot([timeMD[0],timeMD[-1]],[CON3FZ_mean,CON3FZ_mean],color="r")
    t = plt.text((timeMD[0]+timeMD[-1])/2,CON3FZ_mean,f"mean = {CON3FZ_mean:5.3f} N")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, timeMD[-1])
    plt.savefig(folder_path+"\\output_CON3FZ_time_history.png") #save as png
    plt.clf()


    # Fairlead tensions

    plt.xlabel('time (s)')
    plt.ylabel('Fairlead 1 T (kN)')
    plt.title('Fairlead Tension 1')
    plt.plot(timeMD,FAIRTEN1)
    t = plt.text((timeMD[0]+timeMD[-1])/2,FAIRTEN1_mean,f"mean = {FAIRTEN1_mean:5.3f} kN")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    s = plt.text((timeMD[0]+timeMD[-1])/2,(FAIRTEN1_mean+FAIRTEN1_max)/2,f"std = {FAIRTEN1_std:5.3f} kN")
    s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    c = plt.text((timeMD[0]+timeMD[-1])/2,FAIRTEN1_max,f"max = {FAIRTEN1_max:5.3f} kN")
    c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, timeMD[-1])
    plt.savefig(folder_path+"\\output_FAIRTEN1_time_history.png") #save as png
    plt.clf()

    plt.xlabel('time (s)')
    plt.ylabel('Fairlead 2 T (kN)')
    plt.title('Fairlead Tension 2')
    plt.plot(timeMD,FAIRTEN2)
    t = plt.text((timeMD[0]+timeMD[-1])/2,FAIRTEN2_mean,f"mean = {FAIRTEN2_mean:5.3f} kN")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    s = plt.text((timeMD[0]+timeMD[-1])/2,(FAIRTEN2_mean+FAIRTEN2_max)/2,f"std = {FAIRTEN2_std:5.3f} kN")
    s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    c = plt.text((timeMD[0]+timeMD[-1])/2,FAIRTEN2_max,f"max = {FAIRTEN2_max:5.3f} kN")
    c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, timeMD[-1])
    plt.savefig(folder_path+"\\output_FAIRTEN2_time_history.png") #save as png
    plt.clf()

    plt.xlabel('time (s)')
    plt.ylabel('Fairlead 3 T (kN)')
    plt.title('Fairlead Tension 3')
    plt.plot(timeMD,FAIRTEN3)
    t = plt.text((timeMD[0]+timeMD[-1])/2,FAIRTEN3_mean,f"mean = {FAIRTEN3_mean:5.3f} kN")
    t.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    s = plt.text((timeMD[0]+timeMD[-1])/2,(FAIRTEN3_mean+FAIRTEN3_max)/2,f"std = {FAIRTEN3_std:5.3f} kN")
    s.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    c = plt.text((timeMD[0]+timeMD[-1])/2,FAIRTEN3_max,f"max = {FAIRTEN3_max:5.3f} kN")
    c.set_bbox(dict(facecolor='w', alpha=0.8, edgecolor='w'))
    plt.xlim(evalTime, timeMD[-1])
    plt.savefig(folder_path+"\\output_FAIRTEN3_time_history.png") #save as png
    plt.clf()
