# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 10:59:22 2021

@author: Giancarlo&Guido

"""

# Import standard math library and input/output handling, system communication, time logging and sub-process executions
import numpy as np
import subprocess
import os
import shutil
import uuid
import time
import math
from datetime import datetime
from pathlib import Path

# Import third-party custom open-source libraries for input/output handling and software execution
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyFAST.input_output.fast_output_file import FASTOutputFile
from mappp_mooring_response import calc_mooring_restoring_matrix
#import QBladeDllInterface.qbladesys as QBlade
from pyQBlade.qblade_input_file import QBladeInputFile
from pyQBlade.qblade_output_file import QBladeOutputFile
from turbclass import TurbModel

# Import custom libraries for platform and moorings modifications and frequency domain analysis
from moormod import moor_config_openfast,moor_config_qblade
from platmod import plat_config_openfast, plat_config_qblade
from timetofreqdomain import freq_domain_data
from postproc_timehistories import plot_func, plot_func_moor


# OpenFAST exe path
FASTexe = os.getcwd() + r'\\openfast_x64.exe'
    
# Creates and modifies the files necessary to simulate a floating wind turbine in time domain and calculates specific performance from simulation
# xx (list) : design variables values, their order MUST correspond to the keys of model['DESVARIABLES']
# model (object of class "TurbModel"): object of the class representing the simulation and the wind turbine model, behaves like a dictionary
# evaltime  (float) : > 0 initial time instant at which the specific performance of the simulation will be evaluated (used to avoid transient effects)
# penaltyValue (float): value returned by this function if constraints are not satisfied (constraints are defined in model['CONSTRAINTS'])
# filepath_template (string) : filepath of the folder in which the model files are present

def eval_Fobj(xx,\
              model,\
              evalTime = 0, penaltyValue = 9999.9,\
              filepath_template = os.getcwd() + '\\sims\\template_input_files'):
    
    f_max = -1
    # Save outputs to print in "Pop_list.txt", first we save the design variables current values
    outputsToPrint = {}
    k = 1
    for i in xx:
      var_name = 'xx_'+str(k)
      outputsToPrint[var_name] = i
      k += 1

    # Start computer time
    eval_init_time = time.time()
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    # Recognize simulation software from "TurbModel" object
    simSoftware = model['SIMSOFTWARE']
    
    # Create new folder for modified files, if the IDFOLDER is 'auto' the function creates a new folder with a random and unique name
    if model['IDFOLDER'] == 'auto':
        id_folder = uuid.uuid4().hex
        mod_folder_name = createModFolder(model['SIMSOFTWARE'],id_folder)
    else:
        id_folder = model['IDFOLDER']
        mod_folder_name = createModFolder(model['SIMSOFTWARE'],id_folder)
    
    # Open first text file to check the correct execution of this function, with starting time, simulation software and folder id
    file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution1'+'.txt', 'a')
    file_object.write("time: %s, recognized simulation software: %s, folder id: %s \n" % (current_time,simSoftware,id_folder))
    file_object.close()
    
    # Create new turbine model (and copy the files needed by the simulation software) in the new folder
    currentTurbModel = createModFiles(model, filepath_template, mod_folder_name)
    
    if currentTurbModel == -1:
        print('Something went wrong in the creation of the new model, the simulation was not performed!')
        return -1

    MD_Data = -1
    sub_data = -1

    # Modify moorings configuration
    if simSoftware == 'OpenFAST':
        MD_Data=moor_config_openfast(xx,
                                     currentTurbModel,
                                     filepath = filepath_template, filepath_mod = mod_folder_name)
        if MD_Data == -1:
            print('Beware, moorings are not present in this simulation')
        
    elif simSoftware == 'QBlade':
        sub_data=moor_config_qblade(xx,
                                    currentTurbModel,
                                    filepath = filepath_template, filepath_mod = mod_folder_name)
        if sub_data == -1:
            print('Beware, moorings are not present in this simulation')
            
    # Modify platform and re-calculate hydrostatics and hydrodynamics, if needed
    capy_init_time = time.time()
    
    if simSoftware == 'OpenFAST':
        HD_data,ED_data,hst_file,K_hst, draft = plat_config_openfast(xx,
                                                                     currentTurbModel, 
                                                                     filepath = filepath_template,
                                                                     filepath_mod = mod_folder_name)   
    elif simSoftware == 'QBlade':
        sub_data = plat_config_qblade(xx, 
                                      currentTurbModel, 
                                      filepath = filepath_template,
                                      filepath_mod = mod_folder_name)
    
    capy_eval_time = time.time()-capy_init_time
    
    file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution2'+'.txt', 'a')
    file_object.write("%s - Capytaine complete in %.1f s \n" % (id_folder,capy_eval_time))
    file_object.close()
    
    # Update paths
    updatePaths(currentTurbModel,mod_folder_name)
        
    # Change initial displacement, if needed ONLY OPENFAST
    if simSoftware == 'OpenFAST':
        if currentTurbModel['OPTIONS']['FixInitDisplacement']:
         ED_data['PtfmSurge'] = currentTurbModel['OPTIONS']['InitDisplacement'][0]
         ED_data['PtfmSway'] = currentTurbModel['OPTIONS']['InitDisplacement'][0]
         ED_data['PtfmHeave'] = currentTurbModel['OPTIONS']['InitDisplacement'][0]
         ED_data['PtfmRoll'] = currentTurbModel['OPTIONS']['InitDisplacement'][0]
         ED_data['PtfmPitch'] = currentTurbModel['OPTIONS']['InitDisplacement'][0]
         ED_data['PtfmYaw'] =  currentTurbModel['OPTIONS']['InitDisplacement'][0]
         ED_data.write(mod_folder_name + '\\' + currentTurbModel['ELSMODFILENAME'])
    
    # Change final time of simulation
    if currentTurbModel['OPTIONS']['TimeDomainSim']:
      if simSoftware == 'OpenFAST':
        fst_data = FASTInputFile(mod_folder_name+"\\"+currentTurbModel['FSTMODFILENAME'])
        fst_data['TMax'] = currentTurbModel['TMAX']
        fst_data.write(mod_folder_name+"\\"+currentTurbModel['FSTMODFILENAME'])
      elif simSoftware == 'QBlade':
        qb_data = QBladeInputFile(mod_folder_name+"\\"+currentTurbModel['SIMMODFILENAME'])
        dt = qb_data['TIMESTEP']
        qb_data['NUMTIMESTEPS'] = round(currentTurbModel['TMAX']/dt)
        qb_data.write(mod_folder_name+"\\"+currentTurbModel['SIMMODFILENAME'])

    # Calculate restoring matrix of mooring system
    K0 = calc_mooring_restoring_matrix(xx,
                                       currentTurbModel,
                                       filepath_template=filepath_template,mapp_template=currentTurbModel['MAPFILENAME'],
                                       filepath_mod=mod_folder_name,mapp_modfile=currentTurbModel['MAPMODFILENAME'])
    
    #print(K0)

    # Check constraints
    
    for i in currentTurbModel['CONSTRAINTS'].keys():
        if i == 'MaxSurgeExcursion':
            # Surge excursion from thrust force
            surge_excursion = currentTurbModel['FIXVARIABLES']['SurgeForce']/K0[0][0]
            
            file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution3'+'.txt', 'a')
            file_object.write("%s - mapp complete, surge: %.3f m \n" % (id_folder,surge_excursion))
            file_object.close()
            
            # If surge exceeds maximum, quit evaluation and return 
            if surge_excursion>currentTurbModel['CONSTRAINTS']['MaxSurgeExcursion']:
                f_max = 1/penaltyValue
                file_object = open(os.getcwd() + '\\sims\\'+'Pop_list'+'.txt', 'a') 
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                file_object.write("%s - surge out of boundary - %.2f %.4f %.2f %s\n" % (current_time, xx[0], xx[1], surge_excursion,id_folder))
                file_object.close()
                return f_max            

    if simSoftware == 'OpenFAST':
        if currentTurbModel['OPTIONS']['TimeDomainSim']:
         # Define simulation file for OpenFAST
         mainfilemod = mod_folder_name+'\\'+currentTurbModel['FSTMODFILENAME']
         new_mainfilemod='"'+mainfilemod+'"'

         # Run openFAST
         openfast_init_time = time.time()
         proc=subprocess.Popen(FASTexe+" "+ new_mainfilemod,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
         stdout,stderr = proc.communicate()

         # Print OpenFAST output
         print(stdout)
         print(stderr)

         openfast_eval_time = time.time()-openfast_init_time

         # Write evaluation
         file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution4'+'.txt', 'a')
         file_object.write("%s - OpenFAST and costs complete in %.1f s \n" % (id_folder,openfast_eval_time))
         file_object.close()

         output_filename = r''+mod_folder_name+'\\'+currentTurbModel['FSTMODFILENAME'][:-4]+'.outb'

         ## Plot output
         file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution5'+'.txt', 'a')
         file_object.write("%s - output filename %s \n" % (id_folder,output_filename))
         file_object.close()
         outdata=FASTOutputFile(output_filename).toDataFrame()         

        if currentTurbModel['OPTIONS']['EvalCosts']:
         # Get Costs to print
         costs = get_costs(xx,currentTurbModel,filepath_template)
         for i in costs.keys():
            outputsToPrint[i] = costs[i]
        
        if currentTurbModel['OPTIONS']['TimeDomainSim']:
         if currentTurbModel['OPTIONS']['FFTAnalysis']:
            if outdata.loc[:,['Time_[s]']].to_numpy()[-1,0] == model['TMAX']:
               
                # Frequency domain analysis
                freqs_Yaw,Pxx1_Yaw,fft1_Yaw,fft1_normalized_Yaw = freq_domain_data(outdata.loc[outdata['Time_[s]'] > evalTime,['PtfmYaw_[deg]']].to_numpy()[:,0],outdata.loc[outdata['Time_[s]'] > evalTime,['Time_[s]']].to_numpy()[:,0], None,folder_path=mod_folder_name,plot_flag=True,plot_name = "Yaw")
                freq_PSDpeak_Yaw = freqs_Yaw[np.argmax(Pxx1_Yaw)]
                FFTpeak_Yaw = max(fft1_normalized_Yaw)

                freqs_Pitch,Pxx1_Pitch,fft1_Pitch,fft1_normalized_Pitch = freq_domain_data(outdata.loc[outdata['Time_[s]'] > evalTime,['PtfmPitch_[deg]']].to_numpy()[:,0],outdata.loc[outdata['Time_[s]'] > evalTime,['Time_[s]']].to_numpy()[:,0], None,folder_path=mod_folder_name,plot_flag=True,plot_name = "Pitch")
                freq_PSDpeak_Pitch = freqs_Pitch[np.argmax(Pxx1_Pitch)]
                FFTpeak_Pitch = max(fft1_normalized_Pitch)

                freqs_Roll,Pxx1_Roll,fft1_Roll,fft1_normalized_Roll = freq_domain_data(outdata.loc[outdata['Time_[s]'] > evalTime,['PtfmRoll_[deg]']].to_numpy()[:,0],outdata.loc[outdata['Time_[s]'] > evalTime,['Time_[s]']].to_numpy()[:,0], None,folder_path=mod_folder_name,plot_flag=True,plot_name = "Roll")
                freq_PSDpeak_Roll = freqs_Roll[np.argmax(Pxx1_Roll)]
                FFTpeak_Roll = max(fft1_normalized_Roll)

                if math.isinf(freq_PSDpeak_Roll):
                    freq_PSDpeak_Roll = 0
                if math.isinf(freq_PSDpeak_Pitch):
                    freq_PSDpeak_Pitch = 0
                if math.isinf(freq_PSDpeak_Yaw):
                    freq_PSDpeak_Yaw = 0

                outputsToPrint['FFTPeakRoll'] = FFTpeak_Roll
                outputsToPrint['FFTPeakPitch'] = FFTpeak_Pitch
                outputsToPrint['FFTPeakYaw'] = FFTpeak_Yaw
                outputsToPrint['freqPeakRoll'] = freq_PSDpeak_Roll
                outputsToPrint['freqPeakPitch'] = freq_PSDpeak_Pitch
                outputsToPrint['freqPeakYaw'] = freq_PSDpeak_Yaw

                # Objective function evaluation
                f_max = FFTpeak_Yaw
            else:
                f_max = 1/penaltyValue
                file_object = open(os.getcwd() + '\\sims\\'+'Pop_list'+'.txt', 'a') 
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                file_object.write("%s - simulation aborted for fatal error - %.2f %.4f %s \n" % (current_time, xx[0], xx[1] ,id_folder))
                file_object.close()
                return f_max
        
    if simSoftware == 'QBlade':
        
        # Define Output Channels
        outputChannels = [b"X_g COG Pos. [m]",b"Y_g COG Pos. [m]",b"Z_g COG Pos. [m]",b"NP Roll X_l [deg]",b"NP Pitch Y_l [deg]",b"NP Yaw Z_l [deg]",b"Thrust [N]"]

        # Define simulation file for QBlade
        simfilemod = mod_folder_name+'\\'+currentTurbModel['SIMMODFILENAME']
        sim_data = QBladeInputFile(simfilemod)
        timeSteps = sim_data["NUMTIMESTEPS"]
        p = str(Path(simfilemod).resolve())
        currentSimFilePath = str.encode(p)
        # filepath_outputFileName = os.getcwd()
        outputFileName = 'outputSim.outq'
        dst_outputFileName = mod_folder_name+"\\"+outputFileName
        # Run QBlade
        qblade_init_time = time.time()
        simulation = QBlade.QBladeSim()
        simulation.runSimulation(0,32,currentSimFilePath,b'final_project.qpr',dst_outputFileName,timeSteps,outputChannels)

        qblade_eval_time = time.time()-qblade_init_time

        # Get Costs to print
        if currentTurbModel['OPTIONS']['EvalCosts']: 
         costs = get_costs(xx,currentTurbModel,filepath_template)
         for i in costs.keys():
            outputsToPrint[i] = costs[i]
        
        # Write evaluation
        file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution4'+'.txt', 'a')
        file_object.write("%s - QBlade and costs complete in %.1f s \n" % (id_folder,qblade_eval_time))
        file_object.close()
        
        ## Plot output
        file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution5'+'.txt', 'a')
        file_object.write("%s - output filename %s \n" % (id_folder,dst_outputFileName))
        file_object.close()
        outdata=QBladeOutputFile(dst_outputFileName).toDataFrame()
        
        if outdata.loc[:,['Time [s]']].to_numpy()[-1,0] == (model['TMAX']-dt):
            # Frequency domain analysis
            freqs_Yaw,Pxx1_Yaw,fft1_Yaw,fft1_normalized_Yaw = freq_domain_data(outdata.loc[outdata['Time [s]'] > evalTime,['NP Yaw Z_l [deg]']].to_numpy()[:,0],outdata.loc[outdata['Time [s]'] > evalTime,['Time [s]']].to_numpy()[:,0], None,folder_path=mod_folder_name,plot_flag=True,plot_name = "Yaw")
            
            freq_PSDpeak_Yaw = freqs_Yaw[np.argmax(Pxx1_Yaw)]
            FFTpeak_Yaw = max(fft1_normalized_Yaw)
            
            freqs_Pitch,Pxx1_Pitch,fft1_Pitch,fft1_normalized_Pitch = freq_domain_data(outdata.loc[outdata['Time [s]'] > evalTime,['NP Pitch Y_l [deg]']].to_numpy()[:,0],outdata.loc[outdata['Time [s]'] > evalTime,['Time [s]']].to_numpy()[:,0], None,folder_path=mod_folder_name,plot_flag=True,plot_name = "Pitch")
            
            freq_PSDpeak_Pitch = freqs_Pitch[np.argmax(Pxx1_Pitch)]
            FFTpeak_Pitch = max(fft1_normalized_Pitch)
    
            freqs_Roll,Pxx1_Roll,fft1_Roll,fft1_normalized_Roll = freq_domain_data(outdata.loc[outdata['Time [s]'] > evalTime,['NP Roll X_l [deg]']].to_numpy()[:,0],outdata.loc[outdata['Time [s]'] > evalTime,['Time [s]']].to_numpy()[:,0], None,folder_path=mod_folder_name,plot_flag=True,plot_name = "Roll")
            
            freq_PSDpeak_Roll = freqs_Roll[np.argmax(Pxx1_Roll)]
            FFTpeak_Roll = max(fft1_normalized_Roll)    
            
            # Function evaluation
            f_max = FFTpeak_Yaw
        else:
            f_max = 1/penaltyValue
            file_object = open(os.getcwd() + '\\sims\\'+'Pop_list'+'.txt', 'a') 
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            file_object.write("%s - simulation aborted for fatal error - %.2f %.4f %s \n" % (current_time, xx[0], xx[1] ,id_folder))
            file_object.close()
            return f_max

    final_eval_time = time.time()-eval_init_time


    # Count succesful executions
    totalDir = 0
    folders_path = os.getcwd() + '\\sims\\'
    count = 0
    
    for base, dirs, files in os.walk(folders_path):
        # print('Searching in : ',base)
        for directories in dirs:
            totalDir += 1
        for x in files:          
            if x.endswith(".outb"):
                count = count + 1               
            elif x.endswith(".outq"):
                count = count + 1

    totalExec = count

    file_object = open(os.getcwd() + '\\sims\\'+'check_parallel_execution6'+'.txt', 'a')
    file_object.write("%s successful execution of individual - number of .outb/.outq files: %d - total evaluation elapsed time %.2f s - f_max: %.4f \n" % (id_folder,totalExec,final_eval_time,f_max))
    file_object.close()


    # Save some important values from time domain simulation
    if currentTurbModel['OPTIONS']['TimeDomainSim']:
      if simSoftware == 'OpenFAST':
          Time=outdata["Time_[s]"]

          PtfmPitch=outdata["PtfmPitch_[deg]"]
          heeling_max=abs(PtfmPitch[Time>Time.iloc[-1]/2].max())
          heeling_mean=abs(PtfmPitch[Time>Time.iloc[-1]/2].mean())
          outputsToPrint['HeelingMax'] = heeling_max
          outputsToPrint['HeelingMean']= heeling_mean

          PtfmSurge=outdata["PtfmSurge_[m]"]
          surge_max=abs(PtfmSurge[Time>Time.iloc[-1]/2].max())
          surge_mean=abs(PtfmSurge[Time>Time.iloc[-1]/2].mean())
          outputsToPrint['SurgeMax'] = surge_max
          outputsToPrint['SurgeMean']= surge_mean

          RotThrust=outdata["RotThrust_[kN]"]
          thrust_max=abs(RotThrust[Time>Time.iloc[-1]/2].max())
          thrust_mean=abs(RotThrust[Time>Time.iloc[-1]/2].mean())
          outputsToPrint['ThrustMax'] = thrust_max
          outputsToPrint['ThrustMean']= thrust_mean

          PtfmYaw=outdata["PtfmYaw_[deg]"]
          yaw_max=PtfmYaw.max()
          outputsToPrint['YawMax']= yaw_max


      elif simSoftware == 'QBlade':
          Time=outdata['Time [s]']
          PtfmPitch=outdata['NP Pitch Y_l [deg]']
          heeling_max=abs(PtfmPitch[Time>Time.iloc[-1]/2].max())
          heeling_mean=abs(PtfmPitch[Time>Time.iloc[-1]/2].mean())
          PtfmYaw=outdata['NP Yaw Z_l [deg]']
          yaw_max=PtfmYaw.max()
          PtfmSurge=outdata['X_g COG Pos. [m]']
          surge_max=abs(PtfmSurge[Time>Time.iloc[-1]/2].max())
          surge_mean=abs(PtfmSurge[Time>Time.iloc[-1]/2].mean())
          RotThrust=outdata["Thrust [N]"]
          thrust_max=abs(RotThrust[Time>Time.iloc[-1]/2].max())/1000
          thrust_mean=abs(RotThrust[Time>Time.iloc[-1]/2].mean())/1000

    # Write output files 'Pop_list'

    file_object = open(os.getcwd() + '\\sims\\'+'Pop_list'+'.txt', 'a') 
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    width = 12
    decimals = 4
    
    file_object.write('{:s} '.format(current_time))
    for value in outputsToPrint.values():
        file_object.write(f"{value:{width}.{decimals}f} ")

    file_object.write('{:s}\n'.format(id_folder))
    file_object.close()
    
    # Plot outputs, if needed
    plot_t_h_flag = False
    plot_t_h_moor_flag = False
    
    if plot_t_h_flag:
        plot_func(simSoftware,evalTime,outdata,folder_path=mod_folder_name)
    if plot_t_h_moor_flag & (simSoftware == 'OpenFAST'):
        output_moor_filename = r''+mod_folder_name+'\\'+currentTurbModel['FSTMODFILENAME'][:-4]+'.MD.out'
        outdata_moor = FASTOutputFile(output_moor_filename).toDataFrame()
        plot_func_moor(simSoftware,evalTime,outdata_moor,folder_path=mod_folder_name)

    # Check heeling constraint, if needed    
    for i in currentTurbModel['CONSTRAINTS'].keys():
        if i == 'MaxHeelAngle':
            if heeling_mean> currentTurbModel['CONSTRAINTS']['MaxHeelAngle']:
                f_max=1/penaltyValue
            
    return f_max


# Helper Functions

def createModFolder(simSoftware,name):
    
    if simSoftware == 'OpenFAST':
        mod_folder_name = os.getcwd() + f"\\sims\\mod_input_files_{name}"
    elif simSoftware == 'QBlade':
        mod_folder_name = os.getcwd() + f"\\sims\\mod_input_files_QBlade_{name}"
            
    try:
        os.mkdir(mod_folder_name)
    except OSError:
        print ("Creation of the directory %s failed" % mod_folder_name)
    else:
        print ("Successfully created the directory %s" % mod_folder_name)
        
    return mod_folder_name

def createModFiles(model,filepath_template,mod_folder_name):
        
    if model['SIMSOFTWARE'] == 'OpenFAST':
        
        try:
            model.addKeyVal('FSTMODFILENAME',model['FSTFILENAME'])
        except:
            print('No "FAST" driver filename provided, OpenFAST simulation will not be performed!')
            return -1
        try:    
            model.addKeyVal('ELSMODFILENAME',model['ELSFILENAME'])
        except:
            print('No "Elastodyn" filename provided...')
        try:
            model.addKeyVal('MRDMODFILENAME',model['MRDFILENAME'])
        except:
            print('No "Moordyn" filename provided...')
        try:
            model.addKeyVal('HYDMODFILENAME',model['HYDFILENAME'])
        except:
            print('No "Hydrodyn" filename provided...')
        try:
            model.addKeyVal('MAPMODFILENAME',model['MAPFILENAME'])
        except:
            print('No "Mapp" filename provided...')
            
            
    elif model['SIMSOFTWARE'] == 'QBlade':

        try:
            model.addKeyVal('SIMMODFILENAME',model['SIMFILENAME'])
        except:
            print('No "sim" filename provided, QBlade simulation will not be performed!')
            return -1
        try:    
            model.addKeyVal('TRBMODFILENAME',model['TRBFILENAME'])
        except:
            print('No turbine main filename provided...')
        try:
            model.addKeyVal('MAINMODFILENAME',model['MAINFILENAME'])
        except:
            print('No substructure main filename provided...')
        try:
            model.addKeyVal('SUBMODFILENAME',model['SUBFILENAME'])
        except:
            print('No substructure filename provided...')
        try:
            model.addKeyVal('MAPMODFILENAME',model['MAPFILENAME'])
        except:
            print('No "Mapp" filename provided...')
       
    if model['SIMSOFTWARE'] == 'OpenFAST':
        
        # Copy .fst template file to new folder
        src_dir=filepath_template+"\\"+model['FSTFILENAME']
        dst_dir=mod_folder_name+"\\"+model['FSTMODFILENAME']
        shutil.copy(src_dir,dst_dir)
        
        try:
            # Copy blades template to mod folders
            src_dir=filepath_template+"\\"+model['BLDFILENAME']
            dst_dir=mod_folder_name+"\\"+model['BLDFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No blades file was found and copied...')
        
        try:    
            # Copy tower template to mod folders
            src_dir=filepath_template+"\\"+model['TWRFILENAME']
            dst_dir=mod_folder_name+"\\"+model['TWRFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No tower file was found and copied')
        
        try:
            # Copy aerodyn template to mod folders
            src_dir=filepath_template+"\\"+model['ARDFILENAME']
            dst_dir=mod_folder_name+"\\"+model['ARDFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No "Aerodyn" file was found and copied')
        
        try:
            # Copy aerodyn template to mod folders
            src_dir=filepath_template+"\\"+model['IFWFILENAME']
            dst_dir=mod_folder_name+"\\"+model['IFWFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No "Inflow Module" file was found and copied')
            
        try:
            # Copy aerodyn folder to mod folders
            src_dir=filepath_template+model['AERFOLDERNAME']
            dst_dir=mod_folder_name+model['AERFOLDERNAME']
            shutil.copytree(src_dir,dst_dir)
        except:
            print('No "Aerodyn" folder was found and copied')
        
        try:
            # Copy servodyn template to mod folders
            src_dir=filepath_template+"\\"+model['SRVFILENAME']
            dst_dir=mod_folder_name+"\\"+model['SRVFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No "Servodyn" file was found and copied')
        
        try:
            # Copy servodyn structural template to mod folders
            src_dir=filepath_template+"\\"+model['STCFILENAME']
            dst_dir=mod_folder_name+"\\"+model['STCFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No "Servodyn" structural file was found and copied')
            
        try:
            # Copy force template to mod folders
            src_dir=filepath_template+"\\"+model['FRCFILENAME']
            dst_dir=mod_folder_name+"\\"+model['FRCFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No prescribed force for "Servodyn" structural file was found and copied')

        try:
            # Copy moment template to mod folders
            src_dir=filepath_template+"\\"+model['MOMFILENAME']
            dst_dir=mod_folder_name+"\\"+model['MOMFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No prescribed moment for "Servodyn" structural file was found and copied')
        
        try:
            # Copy servodyn template to mod folders
            src_dir=filepath_template+"\\"+model['SRVFOLDERNAME']
            dst_dir=mod_folder_name+"\\"+model['SRVFOLDERNAME']
            shutil.copytree(src_dir,dst_dir)
        except:
            print('No "Servodyn" folder was found and copied')
            
    elif model['SIMSOFTWARE'] == 'QBlade':
        
        # Copy FST template to mod folders
        src_dir=filepath_template+"\\"+model['SIMFILENAME']
        dst_dir=mod_folder_name+"\\"+model['SIMMODFILENAME']
        shutil.copy(src_dir,dst_dir)
        
        try:
            # Copy blades template to mod folders
            src_dir=filepath_template+"\\"+model['BLDFILENAME']
            dst_dir=mod_folder_name+"\\"+model['BLDFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No blades file was found and copied...')
        
        try:    
            # Copy tower template to mod folders
            src_dir=filepath_template+"\\"+model['TWRFILENAME']
            dst_dir=mod_folder_name+"\\"+model['TWRFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No tower file was found and copied')
        
        try:
            # Copy turbine template to mod folders
            src_dir=filepath_template+"\\"+model['TRBFILENAME']
            dst_dir=mod_folder_name+"\\"+model['TRBMODFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No turbine file was found and copied')
            
        try:
            # Copy main template to mod folders
            src_dir=filepath_template+"\\"+model['MAINFILENAME']
            dst_dir=mod_folder_name+"\\"+model['MAINMODFILENAME']
            shutil.copy(src_dir,dst_dir)
        except:
            print('No main structure file was found and copied')   
            
        try:
            # Copy aerodyn folder to mod folders
            src_dir=filepath_template+model['AERFOLDERNAME']
            dst_dir=mod_folder_name+model['AERFOLDERNAME']
            shutil.copytree(src_dir,dst_dir)
        except:
            print('No aerodyn folder was found and copied')
    
    newModel = model
    
    return newModel

def updatePaths(model,mod_folder_name):
    
    if model['SIMSOFTWARE'] == 'OpenFAST':

        fst_data = FASTInputFile(mod_folder_name+"\\"+model['FSTMODFILENAME'])
        try:
            fst_data['EDFile'] = model['ELSMODFILENAME']
        except:
            print('No "Elastodyn" modified file, going to deactivate "Elastodyn" module in OpenFAST simulation...')
            fst_data['CompElast'] = 0
        try:
            fst_data['MooringFile'] = model['MRDMODFILENAME']
        except:
            print('No "Moordyn" modified file, going to deactivate "Moordyn" module in OpenFAST simulation...')
            fst_data['CompMooring'] = 0
        try:
            fst_data['HydroFile'] = model['HYDMODFILENAME']
        except:
            print('No "Hydrodyn" modified file, going to deactivate "Hydrodyn" module in OpenFAST simulation...')
            fst_data['CompHydro'] = 0
        try:
            fst_data['AeroFile'] = model['ARDFILENAME']
        except:
            print('No "Aerodyn" file, going to deactivate "Aerodyn" module in OpenFAST simulation...')
            fst_data['CompAero'] = 0
        try:
            fst_data['InflowFile'] = model['IFWFILENAME']
        except:
            print('No "Inflow Module" file, going to deactivate "Inflow Module" in OpenFAST simulation...')
            fst_data['CompInflow'] = 0
        
        fst_data.write(mod_folder_name+"\\"+model['FSTMODFILENAME'])
        
    elif model['SIMSOFTWARE'] =='QBlade':

        sim_data = QBladeInputFile(mod_folder_name+"\\"+model['SIMMODFILENAME'])
        trb_data = QBladeInputFile(mod_folder_name+"\\"+model['TRBMODFILENAME'])
        main_data = QBladeInputFile(mod_folder_name+"\\"+model['MAINMODFILENAME'])
        sim_data['TURBFILE'] = model['TRBMODFILENAME']
        trb_data['STRUCTURALFILE'] = model['MAINMODFILENAME']
        main_data['SUBFILE'] = model['SUBMODFILENAME']
        sim_data.write(mod_folder_name+"\\"+model['SIMMODFILENAME'])
        trb_data.write(mod_folder_name+"\\"+model['TRBMODFILENAME'])
        main_data.write(mod_folder_name+"\\"+model['MAINMODFILENAME'])
        
    return

# Custom function for calculating Costs of Chains and Horizontal Legs
def get_costs(xx,model,filepath):

    costs = {}                
    MBL = 18700 # kN
    steel_cost = 5 # euro/kg
    horleg_mass = 637.686*1000 # kg weight of 3 braces
    horleg_thick = 0.0564 # m
    rho_steel = 7850 # kg/m^3
    LineNumber = 3
    SparNumber = 3

    a = 0
    
    for i in model['DESVARIABLES'].keys():
        if i == 'AnchorRadius':
            AnchorRadius = model['DESVARIABLES']['AnchorRadius']
            PresAnchorRadius = xx[a]
            a+=1
        elif i == 'LineLengthFactor':
            LineLengthFactor = model['DESVARIABLES']['LineLengthFactor']
            PresLineLengthFactor = xx[a]
            a+=1
        elif i == 'SparDistance':
            SparDistance = model['DESVARIABLES']['SparDistance']
            PresSparDistance = xx[a]
            a+=1
        elif i == 'FairleadRadius':
            FairleadRadius = model['DESVARIABLES']['FairleadRadius']
            PresFairleadRadius = xx[a]
            a+=1
        elif i == 'FairleadHeight':
            FairleadHeight = model['DESVARIABLES']['FairleadHeight']
            PresFairleadHeight = xx[a]
            a+=1
        else:
            a+=1
    
    if 'AnchorRadius' not in model['DESVARIABLES'].keys():
        AnchorRadius = model['FIXVARIABLES']['AnchorRadius']
        PresAnchorRadius = model['FIXVARIABLES']['AnchorRadius']
    if 'FairleadRadius' not in model['DESVARIABLES'].keys():
        FairleadRadius = model['FIXVARIABLES']['FairleadRadius']
        PresFairleadRadius = model['FIXVARIABLES']['FairleadRadius']
    if 'FairleadHeight' not in model['DESVARIABLES'].keys():
        FairleadHeight = model['FIXVARIABLES']['FairleadHeight']
        PresFairleadHeight = model['FIXVARIABLES']['FairleadHeight']
    if 'LineLengthFactor' not in model['DESVARIABLES'].keys():
        LineLengthFactor = model['FIXVARIABLES']['LineLengthFactor']
        PresLineLengthFactor = model['FIXVARIABLES']['LineLengthFactor']

    
    if 'SparDistance' not in model['DESVARIABLES'].keys():
       if 'SparDistance' not in model['FIXVARIABLES'].keys():
           SparDistance = 0
           PresSparDistance = 0
       else:    
           SparDistance = model['FIXVARIABLES']['SparDistance']
           PresSparDistance = model['FIXVARIABLES']['SparDistance']
    
    SparNumber = 3 # number of spars
    
    # To get "Depth" directly from OpenFast or QBlade files
    if model['SIMSOFTWARE'] == 'OpenFAST':
        HD_file=filepath +"\\"+ model['HYDFILENAME']
        HD_Data = FASTInputFile(HD_file)
        WtrDpth=HD_Data['WtrDpth']
    elif model['SIMSOFTWARE'] == 'QBlade':
        sub_file = filepath +"\\"+ model['SUBFILENAME']
        sub_Data = QBladeInputFile(sub_file)
        WtrDpth=sub_Data['WATERDEPTH']

    if 'MoorCosts' in model['OPTIONS']['Costs']:
      # Chain Costs - Chain cost function from "Castillo 2020" who cited C. Consortium
      ref_fairlead_anchor_distance = np.sqrt((WtrDpth+FairleadHeight)**2+(AnchorRadius-FairleadRadius)**2)
      fairlead_anchor_distance = np.sqrt((WtrDpth+PresFairleadHeight)**2+(PresAnchorRadius-PresFairleadRadius)**2)
      ref_chain_cost = LineNumber*(0.0591*MBL-89.69)*ref_fairlead_anchor_distance*LineLengthFactor
      pres_chain_cost = LineNumber*((0.0591*MBL-89.69)*fairlead_anchor_distance*PresLineLengthFactor)
      costs['RefChainCost']=ref_chain_cost
      costs['CurrentChainCost']=pres_chain_cost
    
    if 'BracesCosts' in model['OPTIONS']['Costs']:
      # Horizontal Legs Costs - From Giancarlo Troise 2022 
      ref_horleg_length = SparDistance
      ref_horleg_cost = steel_cost * (horleg_mass + SparNumber*4*np.sqrt(ref_horleg_length/1.22)*horleg_thick*(ref_horleg_length-ref_horleg_length)*rho_steel)
      pres_horleg_cost = steel_cost * (horleg_mass + SparNumber*4*np.sqrt(PresSparDistance/1.22)*horleg_thick*(PresSparDistance-ref_horleg_length)*rho_steel)
      costs['RefBracesCost']=ref_horleg_cost
      costs['CurrentBracesCost']=pres_horleg_cost    
    return costs

if __name__ == '__main__':
   
#   # SOFTWIND CASE
#   templateModel = TurbModel(r'.\sims\template_input_files_QBlade\_DTU10MWSoftwind_modeldefinition_QB.dat')

#   # Possible Types: "TripleSpar", "Spar"
#   templateModel.addKeyVal('PLATFORMTYPE','Spar') 
#   templateFolder = r'.\sims\template_input_files_QBlade'
#
#   # Possible key values for "Spar" model:
#   #        FairleadRadius, FairleadHeight, LineLengthFactor, AnchorRadius (both fixed and design)
#   #        SurgeForce, Depth, LineNumber (fixed only)
#   templateModel.addKeyVal('DESVARIABLES',{'FairleadRadius' : 18, 'FairleadHeight': -15 , 'LineLengthFactor' : 1.05})
#   templateModel.addKeyVal('FIXVARIABLES',{'LineNumber' : 3, 'AnchorRadius': 638.4,'SurgeForce': 1650000, 'Depth' : 200})
#   
#   # Surge constraint is checked prior to the simulation (MAP++); heel angle is checked after the simulation, max between 
#   templateModel.addKeyVal('CONSTRAINTS',{'MaxSurgeExcursion' : 50.0 , 'MaxHeelAngle' : 10}) 
#   
#   templateModel.addKeyVal('TMAX',180) # final simulation time
#   templateModel.addKeyVal('IDFOLDER','auto') # automatic folder name 
#   templateModel.addKeyVal('SIMSOFTWARE','QBlade') # software choice ("OpenFAST" or "QBlade")
#
#   xx = [27,-13.5044,1.07] # single run design variable assignment
#
#   f_max_try = eval_Fobj(xx,templateModel,evalTime = 0,filepath_template = templateFolder)
    
  ## TRIPLE SPAR CASE ##
  # Simulation options:
  templateModel = TurbModel(r'.\sims\template_input_files\_DTU10MW3Spar_modeldefinition.dat')
  templateFolder = r'.\sims\template_input_files'
  templateModel.addKeyVal('OPTIONS',{'TimeDomainSim':True,'FFTAnalysis':True,'EvalCosts':True,'Costs':['MoorCosts','BracesCosts'],\
                                     'FixInitDisplacement':False,'InitDisplacement':[0,0,0,0,0,0]})    
  templateModel.addKeyVal('DESVARIABLES',{'SparDistance': 32.0,'LineLengthFactor' : 1.05676})
  templateModel.addKeyVal('FIXVARIABLES',{'LineNumber' : 3,
                                          'AnchorRadius' : 600.0,
                                          'FairleadRadius' : 54.48, 'FairleadHeight': 8.7,'FairleadDistance': 7.5 + 21.0,
                                          'SparRadius' : 7.5, 'SparHeight' : 53.964, 'Depth':180,
                                          'HPRadius' : 11.25, 'HPHeight' : 0.5,
                                          'SurgeForce': 1650000})
  templateModel.addKeyVal('CONSTRAINTS',{})
  #templateModel.addKeyVal('CONSTRAINTS',{'MaxSurgeExcursion' : 50.0 , 'MaxHeelAngle' : 5})
  templateModel.addKeyVal('PLATFORMTYPE','TripleSpar')
  templateModel.addKeyVal('TMAX',600)
  templateModel.addKeyVal('IDFOLDER','PROVANUOVOPYFAST')
  templateModel.addKeyVal('SIMSOFTWARE','OpenFAST')
  xx = [26.0,1.05676]
  f_max = eval_Fobj(xx,templateModel,evalTime = 0,filepath_template = templateFolder)