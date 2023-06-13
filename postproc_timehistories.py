# -*- coding: utf-8 -*-
"""
Created on Thu May 26 12:06:53 2022

@author: Guido
"""
# Import standard libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pyFAST.input_output.fast_output_file import FASTOutputFile
# Plot time histories 
def plot_func(simSoftware,evalTime,outdata,folder_path='unused'):
    
    if simSoftware == 'OpenFAST':
        time=outdata["Time_[s]"].to_numpy()[:]
        pitch=outdata["PtfmPitch_[deg]"].to_numpy()[:]
        yaw=outdata["PtfmYaw_[deg]"].to_numpy()[:]
        surge=outdata["PtfmSurge_[m]"].to_numpy()[:]
        heave=outdata["PtfmHeave_[m]"].to_numpy()[:]
        roll=outdata["PtfmRoll_[deg]"].to_numpy()[:]
        thrust = outdata["RotThrust_[kN]"].to_numpy()[:]
        power = outdata["RotPwr_[kW]"].to_numpy()[:]
        rotspeed = outdata["RotSpeed_[rpm]"].to_numpy()[:]
        genspeed = outdata["GenSpeed_[rpm]"].to_numpy()[:]
        bladepitch1 = outdata["PtchPMzc1_[deg]"].to_numpy()[:]
        bladepitch2 = outdata["PtchPMzc2_[deg]"].to_numpy()[:]
        bladepitch3 = outdata["PtchPMzc3_[deg]"].to_numpy()[:]

    elif simSoftware == 'QBlade':
        time=outdata['Time [s]'].to_numpy()[:]
        pitch=outdata['NP Pitch Y_l [deg]'].to_numpy()[:]
        yaw=outdata['NP Yaw Z_l [deg]'].to_numpy()[:]
        surge=outdata['X_g COG Pos. [m]'].to_numpy()[:]
        heave=outdata['Z_g COG Pos. [m]'].to_numpy()[:]
        roll=outdata['NP Roll X_l [deg]'].to_numpy()[:]
        
    plt.xlabel('time (s)')
    plt.ylabel('yaw (deg)')
    plt.title('Yaw Time History')
    plt.plot(time,yaw)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_yaw_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('yaw (deg)')
    plt.title('Yaw Complete Time History')
    plt.plot(time,yaw)
    plt.xlim(0, time[-1])
    plt.savefig(folder_path+"\\output_yaw_complete_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('roll (deg)')
    plt.title('Roll Time History')
    plt.plot(time,roll)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_roll_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('surge (m)')
    plt.title('Surge Time History')
    plt.plot(time,surge)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_surge_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('pitch (deg)')
    plt.title('Pitch Time History')
    plt.plot(time,pitch)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_pitch_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('heave (m)')
    plt.title('Heave Time History')
    plt.plot(time,heave)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_heave_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('thrust (kN)')
    plt.title('Thrust Time History')
    plt.plot(time,thrust)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_thrust_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('power (kW)')
    plt.title('Power Time History')
    plt.plot(time,power)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_power_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('b. pitch (deg)')
    plt.title('Blade Pitch Time History')
    plt.plot(time,bladepitch1)
    plt.plot(time,bladepitch2)
    plt.plot(time,bladepitch3)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_bladepitch_time_history.png") #save as png
    plt.clf()
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('rot. speed (RPM)')
    plt.title('Rot. Speed Time History')
    plt.plot(time,rotspeed)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_rotspeed_time_history.png") #save as png
    plt.clf()

def plot_func_moor(simSoftware,evalTime,outdata,folder_path='unused'):
    
    if simSoftware == 'OpenFAST':
        time=outdata["Time_[s]"].to_numpy()[:]
        Con1FZ=outdata["CON1FZ_[N]"].to_numpy()[:]/1000
        Con2FZ=outdata["CON2FZ_[N]"].to_numpy()[:]/1000
        Con3FZ=outdata["CON3FZ_[N]"].to_numpy()[:]/1000
        L1N1PZ=outdata["L1N1PZ_[m]"].to_numpy()[:]
        L2N1PZ=outdata["L2N1PZ_[m]"].to_numpy()[:]
        L3N1PZ=outdata["L3N1PZ_[m]"].to_numpy()[:]
        

    plt.xlabel('time (s)')
    plt.ylabel('L1N1PZ (m)')
    plt.title('L1N1PZ Time History')
    plt.plot(time,L1N1PZ)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_L1N1PZ_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('L2N1PZ (m)')
    plt.title('L2N1PZ Time History')
    plt.plot(time,L2N1PZ)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_L2N1PZ_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('L3N1PZ (m)')
    plt.title('L3N1PZ Time History')
    plt.plot(time,L3N1PZ)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_L3N1PZ_time_history.png") #save as png
    plt.clf()
    
    plt.xlabel('time (s)')
    plt.ylabel('Anch. Fz (kN)')
    plt.title('Con1FZ Time History')
    plt.plot(time,Con1FZ)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_Con1FZ_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('Anch. Fz (kN)')
    plt.title('Con2FZ Time History')
    plt.plot(time,Con2FZ)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_Con2FZ_time_history.png") #save as png
    plt.clf()
    plt.xlabel('time (s)')
    plt.ylabel('Anch. Fz (kN)')
    plt.title('Con3FZ Time History')
    plt.plot(time,Con3FZ)
    plt.xlim(evalTime, time[-1])
    plt.savefig(folder_path+"\\output_Con3FZ_time_history.png") #save as png
    plt.clf()

if __name__ == '__main__':

    output_filename= r'C:\Users\Utente\Desktop\Work_Guido\Projects\FLOATECH\WP4\triple_spar_optimization_helix\mod_input_files_triple_spar_optimized_constforce_helix\DTU10MW3Spar_Param_MPHC_par.outb'
    #output_filename2 = r'C:\Users\Utente\Desktop\Work_Guido\Projects\FLOATECH\WP4\softwind_optimization_helix_0404\optimization_2\mod_input_files_f7b43d5b98a64438bc2df2a460a926bb\SW_SPAR_STATIC_EQUILIBRIUM_EQ3MOORING.MD.out'
    
    templateFolder = r'C:\Users\Utente\Desktop\Work_Guido\Projects\FLOATECH\WP4\triple_spar_optimization_helix\mod_input_files_triple_spar_optimized_constforce_helix'
    simSoftware = 'OpenFAST'
    outdata=FASTOutputFile(output_filename).toDataFrame()
    #outdata_moor=FASTOutputFile(output_filename2).toDataFrame()
    evalTime = 0
    plot_func(simSoftware,evalTime,outdata,templateFolder)
    #plot_func_moor(simSoftware,evalTime,outdata_moor,templateFolder)