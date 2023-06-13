# -*- coding: utf-8 -*-

"""
Created on Mon Oct 25 10:58:10 2021
Optimize floating platform of FOWTs for specific objectives 
@author: Giancarlo Troise&Guido Lazzerini
"""

# Import standard libraries
import scipy as sp
import numpy as np
import time
import os
from datetime import datetime

# Import third-party software input/output software handling
from pyFAST.input_output.fast_input_file import FASTInputFile
from pyQBlade.qblade_input_file import QBladeInputFile
from turbclass import TurbModel

# Import user defined simulation of FOWT specific library
import simFOWT

# Define template model file (containing all subfile names) and folder

templateModel = TurbModel(r'.\sims\template_input_files\_DTU10MWSoftwind_modeldefinition.dat')
templateFolder = r'.\sims\template_input_files'


templateModel.addKeyVal('DESVARIABLES',{'FairleadRadius' : 20, 'FairleadHeight': -20 , 'LineLengthFactor' : 1.05})
templateModel.addKeyVal('FIXVARIABLES',{'LineNumber' : 3, 'AnchorRadius': 638.4,'SurgeForce': 1650000, 'Depth' : 200})
templateModel.addKeyVal('CONSTRAINTS',{'MaxSurgeExcursion' : 100.0 , 'MaxHeelAngle' : 10})
templateModel.addKeyVal('PLATFORMTYPE','Spar')

templateModel.addKeyVal('BOUNDARIES',[(templateModel['DESVARIABLES']['FairleadRadius']*0.75, templateModel['DESVARIABLES']['FairleadRadius']*1.25),
                                      (templateModel['DESVARIABLES']['FairleadHeight']*0.975, templateModel['DESVARIABLES']['FairleadHeight']*1.025),
                                      (templateModel['DESVARIABLES']['LineLengthFactor']*0.975, templateModel['DESVARIABLES']['LineLengthFactor']*1.025)])

templateModel.addKeyVal('TMAX',1200)
templateModel.addKeyVal('IDFOLDER','auto')

# Simulation parameters
evalTime = 600 # simulation starting evaluation time - [s]
penaltyValue = 9999.9 # penalty value for the objective function [-]

#--OLD-- Variables
#LineNumber = 3 # number of mooring lines - [-]
#FairleadRadius = 54.48 # fairlead to Z axis distance - [m]
#FairleadHeight = 8.7 # fairlead to X-Y plane distance - [m]
#FairleadDistance = 7.5 + 21.0 # distance of fairlead from spars axis - [m]
#externalSurgeForce = 1650000.0  # thrust force - [N]

#--OLD-- Constraints
#max_surge_excursion = 25.0 # usually 15% of depth - [m]
#max_heeling_angle = 5  # always around 5° to 10° - [deg]

# Target function definition
def f_target(xx):

    time.sleep(0.25)
    
    # Weights of the optimization
    w_freq = 0.90
    w_cost = 0.10
        
    # Get costs of chains and horizontal legs
    costs = simFOWT.get_costs(xx,templateModel,templateFolder)
    
    # Get yaw amplitude in forced oscillations
    yaw_amp = simFOWT.eval_Fobj(xx,
                                 templateModel,
                                 evalTime,penaltyValue,\
                                 filepath_template = templateFolder)
    
    fval = w_freq*np.abs((1/yaw_amp))\
         + w_cost*(np.abs(costs[2]-costs[0])/costs[0]+np.abs(costs[3]-costs[1])/costs[1])
        
    return fval

# Function to check differential evolution progress at each iteration
def my_callback(xk, convergence, f_val):
    file_object = open('DE_output.txt', 'a') 
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    file_object.write("%s %f %f %f %f %f \n" % (current_time, xk[0], xk[1], xk[2], convergence, f_val))
    file_object.close()

# Optimization process
def main():
        
    # Common operations for each trial configuration of DE #
    t0 = time.perf_counter()

    # Define boundaries for design variables
    boundaries = templateModel['BOUNDARIES']
    x0 = list(templateModel['DESVARIABLES'].values())

    # Algorithm parameters
    n_workers = 2
    n_iters = 2
    n_popsize = 2
        
    # Launch optimization
    res=sp.optimize.differential_evolution(f_target, boundaries, args=(), \
                                           strategy='best1bin', maxiter=n_iters, \
                                           popsize=n_popsize, tol=0.0001, mutation=(0.8, 1.3), \
                                           recombination=0.75, seed=None, callback=my_callback, disp=True, \
                                           polish=False, init='latinhypercube', atol=0, updating='deferred', \
                                           workers=n_workers, \
                                           constraints=(), \
                                           x0=x0)
    
    print('Optimization finished succesfully')

    # Final post-processing
    t_fin = time.perf_counter() - t0
    
    totalDir = 0
    folders_path = os.getcwd() + '\\sims\\'
    count_exec = 0
    count_eval = -1

    for base, dirs, files in os.walk(folders_path):
        # print('Searching in : ',base)
        for directories in dirs:
            totalDir += 1
        for x in files:    
            if x.endswith(".outb"):
                count_exec = count_exec + 1
            elif x.endswith(".outq"):
                count_exec = count_exec + 1
                
    for base, dirs, files in os.walk(folders_path):
        # print('Searching in : ',base)
        for directories in dirs:
            totalDir += 1
        for x in files:    
            if x.endswith(".fst"):
                count_eval = count_eval + 1
            elif x.endswith(".sim"):
                count_eval = count_eval + 1
               
    avg_eval_time = t_fin/count_eval
    
    now = datetime.now()
    current_time = now.strftime("%H.%M.%S")
    
    file_object = open('check_opt_complete_@'+current_time+'.txt', 'a')
    file_object.write("--- Optimization complete @ %s --- \nworkers: %d \nnumber of iterations: %d \npopulation: %d \ntotal evaluations of obj. function: %d \ntotal executions of OpenFAST/QBlade: %d \ntotal elapsed time: %.2f s\naverage evaluation time: %.2f s" % (current_time,n_workers,n_iters,n_popsize, count_eval,count_exec,t_fin,avg_eval_time))
    file_object.close()
                              
    return res

if __name__ == '__main__':
    res=main()
    