# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 18:00:47 2022

@author: Giancarlo
"""

import map_plus_plus.python_driver.mapsys as MAPpp
import numpy as np

def calc_mooring_restoring_action(LineNumber=3, FairleadRadius=54, FairleadHeight=8.7, AnchorRadius=600, LineLengthFactor=1.05, Depth = 200,\
                                 surge=0.0, sway=0.0, heave=0.0, roll=0.0, pitch=0.0 ,yaw=0.0, \
                                 print_plot=False, filepath_template = "unused",mapp_template="unused", \
                                 filepath_mod="unused", mapp_modfile="unused"):
    
    depth = Depth

    # setup MAP 
    mooring_1 = MAPpp.Map()
    
    mooring_1.map_set_sea_depth(depth)      # m
    mooring_1.map_set_gravity(9.81)       # m/s^2
    mooring_1.map_set_sea_density(1025.0) # kg/m^3
        
    # read and modify template for MAPpp input file
    mappp_inputfile = open(filepath_template+"\\"+mapp_template, "r")
    mpp_lines = mappp_inputfile.readlines()
    mappp_inputfile.close()
    
    xf=np.zeros((3)); yf=np.zeros((3)); zf=np.zeros((3))
    xv=np.zeros((3)); yv=np.zeros((3)); zv=np.zeros((3))

    # fixed nodes (anchors)
    xf[0]= AnchorRadius;                           yf[0]=0.0;                                      zf[0]=-depth  # z equal to depth
    xf[1]= AnchorRadius*np.cos(np.deg2rad(120));   yf[1]= AnchorRadius*np.sin(np.deg2rad(120));    zf[1]=-depth  # z equal to depth
    xf[2]= AnchorRadius*np.cos(np.deg2rad(120));   yf[2]=-AnchorRadius*np.sin(np.deg2rad(120));    zf[2]=-depth  # z equal to depth
    #vessel nodes (fairleads)
    xv[0]= FairleadRadius;                         yv[0]=0.0;                                      zv[0]=FairleadHeight  
    xv[1]= FairleadRadius*np.cos(np.deg2rad(120)); yv[1]= FairleadRadius*np.sin(np.deg2rad(120));  zv[1]=FairleadHeight
    xv[2]= FairleadRadius*np.cos(np.deg2rad(120)); yv[2]=-FairleadRadius*np.sin(np.deg2rad(120));  zv[2]=FairleadHeight
    
    mappp_inputfile = open(filepath_mod+"\\"+mapp_modfile, "w")
    #             Node   Type       X        Y         Z            M      B      FX    FY    FZ
    #             (-)    (-)       (m)      (m)       (m)          (kg)   (mˆ3)   (N)   (N)   (N)
    mpp_lines[7] =("1	fix    {:>12.5f} {:>12.5f}  {:>12.5f}       0	    0	   #     #     #	\n".format(xf[0], yf[0],zf[0])).lstrip()
    mpp_lines[8] =("2	fix    {:>12.5f} {:>12.5f}  {:>12.5f}       0	    0	   #     #     #	\n".format(xf[1], yf[1],zf[1])).lstrip()
    mpp_lines[9] =("3	fix    {:>12.5f} {:>12.5f}  {:>12.5f}       0	    0	   #     #     #	\n".format(xf[2], yf[2],zf[2])).lstrip()
    mpp_lines[10]=("4	vessel {:>12.5f} {:>12.5f}  {:>12.5f}       0	    0	   #     #     #	\n".format(xv[0], yv[0], zv[0])).lstrip()
    mpp_lines[11]=("5	vessel {:>12.5f} {:>12.5f}  {:>12.5f}       0	    0	   #     #     #	\n".format(xv[1], yv[1], zv[1])).lstrip()
    mpp_lines[12]=("6	vessel {:>12.5f} {:>12.5f}  {:>12.5f}       0	    0	   #     #     #	\n".format(xv[2], yv[2], zv[2])).lstrip()
    
    # line lengths
    L1=LineLengthFactor*np.sqrt((xf[0]-xv[0])**2+(yf[0]-yv[0])**2+(zf[0]-zv[0])**2)
    L2=LineLengthFactor*np.sqrt((xf[1]-xv[1])**2+(yf[1]-yv[1])**2+(zf[1]-zv[1])**2)
    L3=LineLengthFactor*np.sqrt((xf[2]-xv[2])**2+(yf[2]-yv[2])**2+(zf[2]-zv[2])**2)
    
    #              Line           LineType         UnstrLen        NodeAnch        NodeFair         Outputs               
    #               (-)             (-)               (m)             (-)             (-)             (-)             (-)             (-)            
    mpp_lines[16]=(" 1	            main	          %.3f	           1	           4	        TENSION_ANCH    tension_fair  x_excursion  Z_EXCURSION  ALTITUDE ALTITUDE_ANCH\n" % (L1)).lstrip()
    mpp_lines[17]=(" 2	            main	          %.3f	           2	           5	        TENSION_ANCH    tension_fair  x_excursion  Z_EXCURSION  ALTITUDE ALTITUDE_ANCH\n" % (L2)).lstrip()
    mpp_lines[18]=(" 3	            main	          %.3f	           3	           6	        TENSION_ANCH    tension_fair  x_excursion  Z_EXCURSION  ALTITUDE ALTITUDE_ANCH\n" % (L3)).lstrip()
    
    for element in mpp_lines:
        mappp_inputfile.write(element)
        
    mappp_inputfile.close()
    
    #### run MAPpp #####
    # mooring_1.read_file("./test_input_map.dat") 
    mooring_1.read_file(filepath_mod+"\\"+mapp_modfile)
    mooring_1.summary_file(filepath_mod+"\\"'MAPpp_summary_file.txt')

    mooring_1.init( )
    
    epsilon = 1e-3 # finite difference epsilon
    
    K0 = mooring_1.linear(epsilon)  # restoring matrix in undisplaced position
    
    #displace platform
    mooring_1.displace_vessel(surge,sway,heave,roll,pitch,yaw)
    mooring_1.update_states(0.0,0)
 
    K = mooring_1.linear(epsilon)  # restoring matrix in displaced position

    # We need to call update states after linearization to find the equilibrium
    mooring_1.update_states(0.0,0)
    
    # get forces
    fairlead_H=np.zeros(3); fairlead_V=np.zeros(3);
    fairlead_fx=np.zeros(3); fairlead_fy=np.zeros(3); fairlead_fz=np.zeros(3)
    surge_reaction = 0.0
    sway_reaction  = 0.0
    heave_reaction = 0.0
    roll_reaction  = 0.0
    pitch_reaction = 0.0
    yaw_reaction   = 0.0
    
    for line_number in range(3):
        fairlead_H[line_number],fairlead_V[line_number] = mooring_1.get_fairlead_force_2d(line_number)    
        fairlead_fx[line_number],fairlead_fy[line_number],fairlead_fz[line_number] = mooring_1.get_fairlead_force_3d(line_number)  
        surge_reaction = surge_reaction + fairlead_fx[line_number]
        sway_reaction  = sway_reaction  + fairlead_fy[line_number]
        heave_reaction = heave_reaction + fairlead_fz[line_number]
    
    if print_plot:
        num_points = 20
        empty_lists = [ [] for _ in range(3) ]
        for i in range(0,mooring_1.size_lines()):
            empty_lists[0].append(mooring_1.plot_x(i, num_points))
            empty_lists[1].append(mooring_1.plot_y(i, num_points))
            empty_lists[2].append(mooring_1.plot_z(i, num_points))
            print(i)
        
        print(empty_lists)
        f = open('coordinates_finalconfig2.txt', 'w')
        for row in empty_lists:
            np.savetxt(f, row)
        f.close()

    mooring_1.end()
    return surge_reaction, sway_reaction, heave_reaction, fairlead_fx, fairlead_fy, fairlead_fz, K0, K

def calc_mooring_restoring_matrix(xx,turbModel, 
                                 surge=0.0, sway=0.0, heave=0.0, roll=0.0, pitch=0.0 ,yaw=0.0, \
                                 print_plot=False,filepath_template="unused",mapp_template="unused",filepath_mod="unused",mapp_modfile="unused"):
    
    a = 0
    
    for i in turbModel['DESVARIABLES'].keys():
        if i == 'AnchorRadius':
            AnchorRadius = xx[a]
            a+=1
        elif i == 'LineLengthFactor':
            LineLengthFactor = xx[a]
            a+=1
        elif i == 'FairleadRadius':
            FairleadRadius = xx[a]
            a+=1
        elif i == 'FairleadHeight':
            FairleadHeight = xx[a]
            a+=1
        else:
            a+=1
    
    if 'AnchorRadius' not in turbModel['DESVARIABLES'].keys():
        AnchorRadius = turbModel['FIXVARIABLES']['AnchorRadius']
    if 'FairleadRadius' not in turbModel['DESVARIABLES'].keys():
        FairleadRadius = turbModel['FIXVARIABLES']['FairleadRadius']
    if 'FairleadHeight' not in turbModel['DESVARIABLES'].keys():
        FairleadHeight = turbModel['FIXVARIABLES']['FairleadHeight']
    if 'LineLengthFactor' not in turbModel['DESVARIABLES'].keys():
        LineLengthFactor = turbModel['FIXVARIABLES']['LineLengthFactor']
            
    LineNumber = turbModel['FIXVARIABLES']['LineNumber']
    Depth = turbModel['FIXVARIABLES']['Depth']

    surge_reaction, sway_reaction, heave_reaction, fairlead_fx, fairlead_fy, fairlead_fz, K0, K = \
    calc_mooring_restoring_action(LineNumber, FairleadRadius, FairleadHeight, AnchorRadius, LineLengthFactor, Depth, \
                                 surge, sway, heave, roll, pitch ,yaw, \
                                 print_plot, filepath_template = filepath_template,mapp_template=mapp_template, \
                                 filepath_mod=filepath_mod, mapp_modfile=mapp_modfile)
    return K0

if __name__ == '__main__': 
    surge_reaction, sway_reaction, heave_reaction, fairlead_fx, fairlead_fy, fairlead_fz, K0, K = \
    calc_mooring_restoring_action(LineNumber=3, FairleadRadius=54-(26-34), FairleadHeight=8.7, AnchorRadius=592.4, LineLengthFactor=1.077, Depth = 200.0, \
                                 surge=0.0, sway=0.0, heave=0.0, roll=0.0, pitch=0.0 ,yaw=0.0, \
                                 print_plot=False)
        