# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 22:57:18 2021

@author: Giancarlo&Guido
"""

# Import standard libraries

import numpy as np
import os
import shutil

# Import third-party library

from pyFAST.input_output.fast_input_file import FASTInputFile
from pyQBlade.qblade_input_file import QBladeInputFile

# Import custom libraries

from preproc_floatplat.floatplatgmshes import create_triple_spar_mesh
from preproc_floatplat.floatplatcapyhydrodyn import create_hydrodyn_database, convert_CAPYtoWAMIT_file
from preproc_floatplat.floatplatmmhydrost import write_hydrost_file


def plat_config_openfast(xx,
                         turbModel,
                         filepath='unused',
                         filepath_mod='unused'):
    
    show_flag = False
    calc_hydro_flag = True

    if turbModel['PLATFORMTYPE']=='TripleSpar':
        
        a = 0
        
        for i in turbModel['DESVARIABLES'].keys():
            if i == 'SparDistance':
                SparDistance = xx[a]
                SparDistance_Reference = turbModel['DESVARIABLES']['SparDistance']
                calc_hydro_flag = True
                a+=1
            else:
                a+=1
        
        if 'SparDistance' not in turbModel['DESVARIABLES'].keys():
         SparDistance = turbModel['FIXVARIABLES']['SparDistance']
         SparDistance_Reference = turbModel['FIXVARIABLES']['SparDistance']
         calc_hydro_flag = False
      
        #read ElastoDyn input template file 
        ED_data = FASTInputFile(filepath+'\\'+turbModel['ELSFILENAME'])
        
        # read Hydrodyn input template file 
        HD_data = FASTInputFile(filepath+'\\'+turbModel['HYDFILENAME'])
      
        hstfile = 0
        K_hst = 0
        draft = 0

        if calc_hydro_flag:

             # Set poisition of spars based on the modified distance for Morison elements
            HD_data["Joints"][0,1]=SparDistance;                    HD_data["Joints"][0,2]=0; 
            HD_data["Joints"][1,1]=SparDistance;                    HD_data["Joints"][1,2]=0; 
            HD_data["Joints"][2,1]=SparDistance;                    HD_data["Joints"][2,2]=0;
            HD_data["Joints"][3,1]=SparDistance*np.cos(2*np.pi/3);  HD_data["Joints"][3,2]=SparDistance*np.sin(2*np.pi/3); 
            HD_data["Joints"][4,1]=SparDistance*np.cos(2*np.pi/3);  HD_data["Joints"][4,2]=SparDistance*np.sin(2*np.pi/3); 
            HD_data["Joints"][5,1]=SparDistance*np.cos(2*np.pi/3);  HD_data["Joints"][5,2]=SparDistance*np.sin(2*np.pi/3);
            HD_data["Joints"][6,1]=SparDistance*np.cos(-2*np.pi/3); HD_data["Joints"][6,2]=SparDistance*np.sin(-2*np.pi/3); 
            HD_data["Joints"][7,1]=SparDistance*np.cos(-2*np.pi/3); HD_data["Joints"][7,2]=SparDistance*np.sin(-2*np.pi/3); 
            HD_data["Joints"][8,1]=SparDistance*np.cos(-2*np.pi/3); HD_data["Joints"][8,2]=SparDistance*np.sin(-2*np.pi/3);

            # Set location of hydrodynamic file (first is the name)
            hydro_folder_name_temp = r'"triple_spar_mesh_mod\triple_spar_mesh_mod"'
            hydro_folder_name = filepath_mod + "\\triple_spar_mesh_mod"
        
            try:
                os.mkdir(hydro_folder_name)
            except OSError:
                print ("Creation of the directory %s failed" % hydro_folder_name)
            else:
                print ("Successfully created the directory %s" % hydro_folder_name)
                
            HD_data["PotFile"] = hydro_folder_name_temp

            spar_radius =  turbModel['FIXVARIABLES']['SparRadius']
            spar_height = turbModel['FIXVARIABLES']['SparHeight']
            hp_radius = turbModel['FIXVARIABLES']['HPRadius']
            hp_thickness = turbModel['FIXVARIABLES']['HPHeight']
            spar1_X = SparDistance*np.cos(0)
            spar1_Y = SparDistance*np.sin(0)
            spar_vertical_divisions = 20
            circle_divisions = 8
            hp_vertical_divisions = 2
            hp_sup_divisions = 4
            hp_inf_divisions = 8
            draft = spar_height + hp_thickness

            input_mesh_file = filepath_mod + "\\"+"triple_spar_mesh_mod.msh"
            create_triple_spar_mesh(input_mesh_file[0:-4],spar_radius,spar_height,hp_radius,hp_thickness,spar1_X,spar1_Y,spar_vertical_divisions, circle_divisions,hp_vertical_divisions,hp_sup_divisions,hp_inf_divisions,show=show_flag)
        
            # omegas in rad/s and direction in rad
            omega_step = 0.10
            omega_start = (2*np.pi)/300
            omega_end = (2*np.pi)/2
        
            omegas = np.arange(omega_start,omega_end,omega_step)
            np.append(omegas, [np.infty], axis=0) # fix the vector for infinite freq
        
            wave_directions = [0]
        
            hd_dataset = create_hydrodyn_database(input_mesh_file, omegas, wave_directions=wave_directions)
            isConverted = convert_CAPYtoWAMIT_file(hd_dataset,hydro_folder_name+"\\triple_spar_mesh_mod")      
        
            # generate hydrostatic file
            output_hst_namefile=hydro_folder_name+"\\triple_spar_mesh_mod"
            is_hstfile_gen, hstfile,K_hst,K_dim=write_hydrost_file(input_mesh_file,output_hst_namefile ,CoG_Z = 0, rho=1023, g=9.81, show=show_flag)

            # change inertia
            PtfmRIner=ED_data["PtfmRIner"]
            PtfmPIner=ED_data["PtfmPIner"]
            PtfmYIner=ED_data["PtfmYIner"]
            PtfmMass=ED_data["PtfmMass"]

            PtfmRIner_mod=PtfmRIner-(2*PtfmMass/3)*(SparDistance_Reference*np.sin(np.pi/3))**2+(2*PtfmMass/3)*(SparDistance*np.sin(np.pi/3))**2
            PtfmPIner_mod=PtfmPIner-(PtfmMass/3)*SparDistance_Reference**2+(PtfmMass/3)*SparDistance**2 \
                                   -(2*PtfmMass/3)*(SparDistance_Reference*np.cos(np.pi/3))**2+(2*PtfmMass/3)*(SparDistance*np.cos(np.pi/3))**2
            PtfmYIner_mod=PtfmYIner-(PtfmMass)*SparDistance_Reference**2+(PtfmMass)*SparDistance**2

            ED_data["PtfmRIner"]=PtfmRIner_mod
            ED_data["PtfmPIner"]=PtfmPIner_mod
            ED_data["PtfmYIner"]=PtfmYIner_mod

        else:
            try:
                # Copy hydrodynamic potential files folder to mod folders
                src_dir=filepath+turbModel['HYDFOLDERNAME']
                dst_dir=filepath_mod+turbModel['HYDFOLDERNAME']
                shutil.copytree(src_dir,dst_dir)
            except:
                print('No hydrodynamic folder was found and copied')
        
        # write hydrodyn and elastodyn modified file
        HD_data.write(filepath_mod+'\\'+turbModel['HYDMODFILENAME'])
        ED_data.write(filepath_mod+'\\'+turbModel['ELSMODFILENAME'])

    elif turbModel['PLATFORMTYPE'] == 'Spar':
              
        #read ElastoDyn input template file 
        ED_data = FASTInputFile(filepath+'\\'+turbModel['ELSFILENAME'])
        
        # read hydrodyn input template file 
        HD_data = FASTInputFile(filepath+'\\'+turbModel['HYDFILENAME'])
            
        try:
            # Copy hydrodynamic potential files folder to mod folders
            src_dir=filepath+turbModel['HYDFOLDERNAME']
            dst_dir=filepath_mod+turbModel['HYDFOLDERNAME']
            shutil.copytree(src_dir,dst_dir)
        except:
            print('No hydrodynamic folder was found and copied')
        
        # write hydrodyn and elastodyn modified file
        HD_data.write(filepath_mod+'\\'+turbModel['HYDMODFILENAME'])
        ED_data.write(filepath_mod+'\\'+turbModel['ELSMODFILENAME'])
        
        hstfile = 0
        K_hst = 0
        draft = 0

    elif turbModel['PLATFORMTYPE'] == 'HydraSpar':
        a = 0
        
        for i in turbModel['DESVARIABLES'].keys():
            if i == 'SparDistance':
                SparDistance = xx[a]
                SparDistance_Reference = turbModel['DESVARIABLES']['SparDistance']
                calc_hydro_flag = True
                a+=1
            else:
                a+=1
        
        if 'SparDistance' not in turbModel['DESVARIABLES'].keys():
         SparDistance = turbModel['FIXVARIABLES']['SparDistance']
         SparDistance_Reference = turbModel['FIXVARIABLES']['SparDistance']
         calc_hydro_flag = False
      
        #read ElastoDyn input template file 
        ED_data = FASTInputFile(filepath+'\\'+turbModel['ELSFILENAME'])
        
        # read Hydrodyn input template file 
        HD_data = FASTInputFile(filepath+'\\'+turbModel['HYDFILENAME'])
      
        hstfile = 0
        K_hst = 0
        draft = 0

        if calc_hydro_flag:

             # Set poisition of spars based on the modified distance for Morison elements
            HD_data["Joints"][0,1]=SparDistance;                    HD_data["Joints"][0,2]=0; 
            HD_data["Joints"][1,1]=SparDistance;                    HD_data["Joints"][1,2]=0; 
            HD_data["Joints"][2,1]=SparDistance;                    HD_data["Joints"][2,2]=0;
            HD_data["Joints"][3,1]=SparDistance*np.cos(2*np.pi/3);  HD_data["Joints"][3,2]=SparDistance*np.sin(2*np.pi/3); 
            HD_data["Joints"][4,1]=SparDistance*np.cos(2*np.pi/3);  HD_data["Joints"][4,2]=SparDistance*np.sin(2*np.pi/3); 
            HD_data["Joints"][5,1]=SparDistance*np.cos(2*np.pi/3);  HD_data["Joints"][5,2]=SparDistance*np.sin(2*np.pi/3);
            HD_data["Joints"][6,1]=SparDistance*np.cos(-2*np.pi/3); HD_data["Joints"][6,2]=SparDistance*np.sin(-2*np.pi/3); 
            HD_data["Joints"][7,1]=SparDistance*np.cos(-2*np.pi/3); HD_data["Joints"][7,2]=SparDistance*np.sin(-2*np.pi/3); 
            HD_data["Joints"][8,1]=SparDistance*np.cos(-2*np.pi/3); HD_data["Joints"][8,2]=SparDistance*np.sin(-2*np.pi/3);

            # Set location of hydrodynamic file (first is the name)
            hydro_folder_name_temp = r'"triple_spar_mesh_mod\triple_spar_mesh_mod"'
            hydro_folder_name = filepath_mod + "\\triple_spar_mesh_mod"
        
            try:
                os.mkdir(hydro_folder_name)
            except OSError:
                print ("Creation of the directory %s failed" % hydro_folder_name)
            else:
                print ("Successfully created the directory %s" % hydro_folder_name)
                
            HD_data["PotFile"] = hydro_folder_name_temp

            spar_radius =  turbModel['FIXVARIABLES']['SparRadius']
            spar_height = turbModel['FIXVARIABLES']['SparHeight']
            hp_radius = turbModel['FIXVARIABLES']['HPRadius']
            hp_thickness = turbModel['FIXVARIABLES']['HPHeight']
            spar1_X = SparDistance*np.cos(0)
            spar1_Y = SparDistance*np.sin(0)
            spar_vertical_divisions = 20
            circle_divisions = 8
            hp_vertical_divisions = 2
            hp_sup_divisions = 4
            hp_inf_divisions = 8
            draft = spar_height + hp_thickness

            input_mesh_file = filepath_mod + "\\"+"triple_spar_mesh_mod.msh"
            create_triple_spar_mesh(input_mesh_file[0:-4],spar_radius,spar_height,hp_radius,hp_thickness,spar1_X,spar1_Y,spar_vertical_divisions, circle_divisions,hp_vertical_divisions,hp_sup_divisions,hp_inf_divisions,show=show_flag)
        
            # omegas in rad/s and direction in rad
            omega_step = 0.10
            omega_start = (2*np.pi)/300
            omega_end = (2*np.pi)/2
        
            omegas = np.arange(omega_start,omega_end,omega_step)
            np.append(omegas, [np.infty], axis=0) # fix the vector for infinite freq
        
            wave_directions = [0]
        
            hd_dataset = create_hydrodyn_database(input_mesh_file, omegas, wave_directions=wave_directions)
            isConverted = convert_CAPYtoWAMIT_file(hd_dataset,hydro_folder_name+"\\triple_spar_mesh_mod")      
        
            # generate hydrostatic file
            output_hst_namefile=hydro_folder_name+"\\triple_spar_mesh_mod"
            is_hstfile_gen, hstfile,K_hst,K_dim=write_hydrost_file(input_mesh_file,output_hst_namefile ,CoG_Z = 0, rho=1023, g=9.81, show=show_flag)

            # change inertia
            PtfmRIner=ED_data["PtfmRIner"]
            PtfmPIner=ED_data["PtfmPIner"]
            PtfmYIner=ED_data["PtfmYIner"]
            PtfmMass=ED_data["PtfmMass"]

            PtfmRIner_mod=PtfmRIner-(2*PtfmMass/3)*(SparDistance_Reference*np.sin(np.pi/3))**2+(2*PtfmMass/3)*(SparDistance*np.sin(np.pi/3))**2
            PtfmPIner_mod=PtfmPIner-(PtfmMass/3)*SparDistance_Reference**2+(PtfmMass/3)*SparDistance**2 \
                                   -(2*PtfmMass/3)*(SparDistance_Reference*np.cos(np.pi/3))**2+(2*PtfmMass/3)*(SparDistance*np.cos(np.pi/3))**2
            PtfmYIner_mod=PtfmYIner-(PtfmMass)*SparDistance_Reference**2+(PtfmMass)*SparDistance**2

            ED_data["PtfmRIner"]=PtfmRIner_mod
            ED_data["PtfmPIner"]=PtfmPIner_mod
            ED_data["PtfmYIner"]=PtfmYIner_mod

        else:
            try:
                # Copy hydrodynamic potential files folder to mod folders
                src_dir=filepath+turbModel['HYDFOLDERNAME']
                dst_dir=filepath_mod+turbModel['HYDFOLDERNAME']
                shutil.copytree(src_dir,dst_dir)
            except:
                print('No hydrodynamic folder was found and copied')
        
        # write hydrodyn and elastodyn modified file
        HD_data.write(filepath_mod+'\\'+turbModel['HYDMODFILENAME'])
        ED_data.write(filepath_mod+'\\'+turbModel['ELSMODFILENAME'])

    return HD_data,ED_data,hstfile,K_hst,draft

def plat_config_qblade(xx,
                       turbModel,
                       filepath='unused',
                       filepath_mod='unused'):
    
    show_flag = False
    calc_hydro_flag = True

    if turbModel['PLATFORMTYPE']=='TripleSpar':
    
        a = 0
        
        for i in turbModel['DESVARIABLES'].keys():
            if i == 'SparDistance':
                SparDistance = xx[a]
                SparDistance_Reference = turbModel['DESVARIABLES']['SparDistance']
                calc_hydro_flag = True
                a+=1
            else:
                SparDistance = turbModel['FIXVARIABLES']['SparDistance']
                SparDistance_Reference = turbModel['FIXVARIABLES']['SparDistance']
                calc_hydro_flag = False
                a+=1
               
        # read substructure input template file
        # right now the file could be already been created by moormod (check if exists), in that case we use the already modified data
        
        if os.path.isfile(filepath_mod+"\\"+turbModel['SUBMODFILENAME']):
            sub_data = QBladeInputFile(filepath_mod+"\\"+turbModel['SUBMODFILENAME'])
        else:
            sub_data = QBladeInputFile(filepath+"\\"+turbModel['SUBFILENAME'])     
              
        #draft = spar_height + hp_thickness
        
        if calc_hydro_flag:

        # Set poisition of spars based on the modified distance for Morison elements
        
            X1 = SparDistance
            X2 = -SparDistance*np.cos(np.pi/3)
            Y2 = SparDistance*np.sin(np.pi/3)
            Y3 = -SparDistance*np.sin(np.pi/3)

            X1_CONN = SparDistance + spar_radius
            X2_CONN = -(SparDistance+ spar_radius)*np.cos(np.pi/3)
            Y2_CONN = (SparDistance+ spar_radius)*np.sin(np.pi/3)
            Y3_CONN = -(SparDistance+ spar_radius)*np.sin(np.pi/3)

            # Change first spar position
            sub_data.modTabElm("SUBJOINTS",range(1,5),"JointX",X1)
            # Change second and third spar positions
            sub_data.modTabElm("SUBJOINTS",range(5,13),"JointX",X2)
            sub_data.modTabElm("SUBJOINTS",range(5,9),"JointY",Y2)
            sub_data.modTabElm("SUBJOINTS",range(9,13),"JointY",Y3)

            # Change first connection to fairlead position
            sub_data.modTabElm("SUBJOINTS",13,"JointX",X1_CONN)
            # Change second and third connections to fairlead position
            sub_data.modTabElm("SUBJOINTS",range(14,16),"JointX",X2_CONN)
            sub_data.modTabElm("SUBJOINTS",14,"JointY",Y2_CONN)
            sub_data.modTabElm("SUBJOINTS",15,"JointY",Y3_CONN)

            # Set location of hydrodynamic file (first is the name)
            
            hydro_folder_name = filepath_mod + r'\HydroData'
            
            relative_path_to_hydro_file_rad = r'\HydroData\triple_spar_mesh_mod.1'
            relative_path_to_hydro_file_exc =  r'\HydroData\triple_spar_mesh_mod.3'
            
            try:
                os.mkdir(hydro_folder_name)
            except OSError:
                print ("Creation of the directory %s failed" % hydro_folder_name)
            else:
                print ("Successfully created the directory %s" % hydro_folder_name)
    
            sub_data["POT_RAD_FILE"] = relative_path_to_hydro_file_rad
            sub_data["POT_EXC_FILE"] = relative_path_to_hydro_file_exc 

            spar_radius =  turbModel['FIXVARIABLES']['SparRadius']
            spar_height = turbModel['FIXVARIABLES']['SparHeight']
            hp_radius = turbModel['FIXVARIABLES']['HPRadius']
            hp_thickness = turbModel['FIXVARIABLES']['HPHeight']
            spar1_X = SparDistance*np.cos(0)
            spar1_Y = SparDistance*np.sin(0)
            spar_vertical_divisions = 20
            circle_divisions = 8
            hp_vertical_divisions = 2
            hp_sup_divisions = 4
            hp_inf_divisions = 8 

            input_mesh_file = filepath_mod + "\\triple_spar_mesh_mod.msh"
            create_triple_spar_mesh(input_mesh_file[0:-4],spar_radius,spar_height,hp_radius,hp_thickness,spar1_X,spar1_Y,spar_vertical_divisions, circle_divisions,hp_vertical_divisions,hp_sup_divisions,hp_inf_divisions,show=show_flag)
        
            # omegas in rad/s and direction in rad
            omega_step = 0.10
            omega_start = (2*np.pi)/300
            omega_end = (2*np.pi)/2
        
            omegas = np.arange(omega_start,omega_end,omega_step)
            np.append(omegas, [np.infty], axis=0) # fix the vector for infinite freq
            wave_directions = [0]
        
            hd_dataset = create_hydrodyn_database(input_mesh_file, omegas, wave_directions)
            isConverted = convert_CAPYtoWAMIT_file(hd_dataset,hydro_folder_name+"\\triple_spar_mesh_mod")
        
            #change delta of wave directions
            sub_data["DELTA_FREQ_RAD"] = str(round(omega_step,2))
            sub_data["DELTA_FREQ_DIFF"] = str(round(omega_step,2))
            #sub_data["DELTA_DIR_DIFF"] = str(round((np.abs(wave_directions[1] - wave_directions[0])),2))

            sub_data["USE_EXCITATION"] = True

            # generate hydrostatic file (USELESS FOR QBlade)
            #output_hst_namefile=hydro_folder_name+"\\triple_spar_mesh_mod"
            #is_hstfile_gen, hstfile,K_hst,K_dim=write_hydrost_file(input_mesh_file,output_hst_namefile ,CoG_Z = 0, rho=1023, g=9.81, show=False)

            # change inertia
            PtfmRIner=sub_data["SUB_MASS"][3][3]
            PtfmPIner=sub_data["SUB_MASS"][4][4]
            PtfmYIner=sub_data["SUB_MASS"][5][5]
            PtfmMass=sub_data["SUB_MASS"][0][0]

            PtfmRIner_mod=PtfmRIner-(2*PtfmMass/3)*(SparDistance_Reference*np.sin(np.pi/3))**2+(2*PtfmMass/3)*(SparDistance*np.sin(np.pi/3))**2
            PtfmPIner_mod=PtfmPIner-(PtfmMass/3)*SparDistance_Reference**2+(PtfmMass/3)*SparDistance**2 \
                                   -(2*PtfmMass/3)*(SparDistance_Reference*np.cos(np.pi/3))**2+(2*PtfmMass/3)*(SparDistance*np.cos(np.pi/3))**2
            PtfmYIner_mod=PtfmYIner-(PtfmMass)*SparDistance_Reference**2+(PtfmMass)*SparDistance**2

            sub_data.modTabElm("SUB_MASS",3,3,PtfmRIner_mod)
            sub_data.modTabElm("SUB_MASS",4,4,PtfmPIner_mod)
            sub_data.modTabElm("SUB_MASS",5,5,PtfmYIner_mod)


        else:
            try:
               # Copy hydrodynamic potential files folder to mod folders
               src_dir=filepath+turbModel['HYDFOLDERNAME']
               dst_dir=filepath_mod+turbModel['HYDFOLDERNAME']
               shutil.copytree(src_dir,dst_dir)
            except:
               print('No hydrodynamic folder was found and copied')

        # write mooring file
        sub_data.write(filepath_mod+"\\"+turbModel['SUBMODFILENAME'])
   
    elif turbModel['PLATFORMTYPE'] == 'Spar':

        if os.path.isfile(filepath_mod+"\\"+turbModel['SUBMODFILENAME']):
            sub_data = QBladeInputFile(filepath_mod+"\\"+turbModel['SUBMODFILENAME'])
        else:
            sub_data = QBladeInputFile(filepath+"\\"+turbModel['SUBFILENAME'])
        
        # Set location of hydrodynamic file (first is the name)
        
        hydro_folder_name = filepath_mod + turbModel['HYDFOLDERNAME']
        
        relative_path_to_hydro_file_rad =  turbModel['HYDFOLDERNAME'] + r'\sw_spar.1'
        relative_path_to_hydro_file_exc =  turbModel['HYDFOLDERNAME'] + r'\sw_spar.3'
        
        try:
            src_dir=filepath+turbModel['HYDFOLDERNAME']
            dst_dir=filepath_mod+turbModel['HYDFOLDERNAME']
            shutil.copytree(src_dir,dst_dir)
        except OSError:
            print ("Creation of the directory %s failed" % hydro_folder_name)
        else:
            print ("Successfully created the directory %s" % hydro_folder_name)

        # write mooring file
        sub_data.write(filepath_mod+"\\"+turbModel['SUBMODFILENAME'])

    return sub_data