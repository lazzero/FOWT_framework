# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 12:34:48 2021

#  Module for processing hydrodynamics for FLOATING PLATFORMS of wind turbines with "Capytaine"
#  Run it as an independent script if you need it.
#  Contains: 
#  - functions: create_hydrodyn_database, convert_CAPYtoWAMIT_file
#
# omegas is in rad/s and wave direction in rad

@author: Guido Lazzerini (99.9%) & Giancarlo Troise (<1%)


"""
import capytaine as cpt
import numpy as np
import matplotlib.pyplot as plt
import os

import xarray as xr

from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO)


def create_hydrodyn_database(input_mesh_file, omegas, wave_directions = [0],*, show=False,depth = np.infty):
    # SOLVE BEM PROBLEMS
    body = cpt.FloatingBody.from_file(input_mesh_file)  # msh file
    body.add_all_rigid_body_dofs()
    body.keep_immersed_part()

    if show : body.show()

    if not np.isin(0, wave_directions):
        print('Error zero direction must be present in wave direction vector')
        return -1
    
    bem_solver = cpt.BEMSolver()
    problems = []
    for wave_direction in wave_directions:
        for omega in omegas:
            if wave_direction==0:
                problems += [cpt.RadiationProblem(sea_bottom=-depth,omega=omega, body=body, radiating_dof=dof) for dof in body.dofs]
            problems += [cpt.DiffractionProblem(sea_bottom=-depth,omega=omega, body=body, wave_direction=wave_direction)]
    results = [bem_solver.solve(problem) for problem in problems]
    #*radiation_results, diffraction_result = results
    dataset = cpt.assemble_dataset(results)

    # save to netCDF format file
    dataset.to_netcdf(input_mesh_file[:-4]+".cpt",engine='h5netcdf',invalid_netcdf=True)

    return dataset

# CONVERT CAPYTAINE OUTPUT TO WAMIT OUTPUT FILES ".1", ".3"

def convert_CAPYtoWAMIT_file(dataset,output_hyd_namefile):
    body_name = dataset.body_name.item()
    body = cpt.FloatingBody.from_file(body_name)  # msh file
    body.add_all_rigid_body_dofs()
    body.keep_immersed_part()
    rho = dataset['rho'].item()
    g = dataset['g'].item()

    omegas = dataset['omega'].values.tolist()
    # fix last value to infinity for 0 period necessary in OpenFast 
    omegas[-1] = np.infty
    wave_directions = dataset['wave_direction'].values.tolist()
    
    # Create name of file
    first_file_ext = '.1'
    second_file_ext = '.3'
    
    first_file = "".join([output_hyd_namefile,first_file_ext])
    second_file = "".join([output_hyd_namefile,second_file_ext])
    
    # datetime object containing current date and time
    now = datetime.now()
    
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
    
    first_file = "".join([output_hyd_namefile,dt_string,first_file_ext]) if os.path.exists(first_file) else first_file
    second_file = "".join([output_hyd_namefile,dt_string,second_file_ext]) if os.path.exists(second_file) else second_file
    
    file_object_1 = open(first_file, "w+")
    file_object_2 = open(second_file, "w+")
    
    # skip header
    for omega_k in omegas:
     for wave_direct_k in wave_directions:
       for dof_i in body.dofs:
        if dof_i == "Surge":
            i = 1
        elif dof_i == "Sway":
            i = 2
        elif dof_i == "Heave":
            i = 3
        elif dof_i == "Roll":
            i = 4
        elif dof_i == "Pitch":
            i = 5
        else:
            i = 6
        if omega_k != np.infty and omega_k != omegas[0]:
            FK_force = dataset['Froude_Krylov_force'].sel(omega = omega_k, wave_direction = wave_direct_k, influenced_dof=dof_i).item()
            D_force = dataset['diffraction_force'].sel(omega = omega_k, wave_direction = wave_direct_k, influenced_dof=dof_i).item()
            adim_exc_force = (FK_force + D_force)/(rho*g)
            mod_aef = np.abs(adim_exc_force)
            pha_aef = np.angle(adim_exc_force)
            rea_aef = np.real(adim_exc_force)
            ima_aef = np.imag(adim_exc_force)
            line_2 = "{period:.3E}   {direct:.3E}   {dof1:n}   {element1:.5E}   {element2:.5E}   {element3:.5E}   {element4:.5E}\n".format(period = (2*np.pi)/omega_k, direct=wave_direct_k*180/np.pi  ,dof1=i, element1=mod_aef, element2=-pha_aef*180/np.pi,element3=rea_aef,element4=-ima_aef)
            file_object_2.write(line_2)
        if wave_direct_k == 0:
         for dof_j in body.dofs:
            if dof_j == "Surge":
                j = 1
            elif dof_j == "Sway":
                j = 2
            elif dof_j == "Heave":
                j = 3
            elif dof_j == "Roll":
                j = 4
            elif dof_j == "Pitch":
                j = 5
            else:
                j = 6
            if omega_k == omegas[0]:
                line_1 = "{period:.3E}   {dof1:n}   {dof2:n}   {element1:.5E}   {element2:.5E}\n".format(period = -1, dof1=i, dof2=j, element1=dataset['added_mass'].sel(omega = omega_k,radiating_dof = dof_i, influenced_dof=dof_j).item()/rho, element2=dataset['radiation_damping'].sel(omega = omega_k,radiating_dof = dof_i, influenced_dof=dof_j).item()/(rho*omega_k))
            elif omega_k == np.infty:
                line_1 = "{period:.3E}   {dof1:n}   {dof2:n}   {element1:.5E}   {element2:.5E}\n".format(period = 0, dof1=i, dof2=j, element1=dataset['added_mass'].sel(omega = dataset.omega.to_numpy()[-1],radiating_dof = dof_i, influenced_dof=dof_j).item()/rho, element2=dataset['radiation_damping'].sel(omega = dataset.omega.to_numpy()[-1],radiating_dof = dof_i, influenced_dof=dof_j).item()/(rho*dataset.omega.to_numpy()[-1]))
            else:
                line_1 = "{period:.3E}   {dof1:n}   {dof2:n}   {element1:.5E}   {element2:.5E}\n".format(period = (2*np.pi)/omega_k, dof1=i, dof2=j, element1=dataset['added_mass'].sel(omega = omega_k,radiating_dof = dof_i, influenced_dof=dof_j).item()/rho, element2=dataset['radiation_damping'].sel(omega = omega_k,radiating_dof = dof_i, influenced_dof=dof_j).item()/(rho*omega_k))
            file_object_1.write(line_1)
            
    file_object_1.close()
    file_object_2.close()
    
    plot_flag = False
    
    if plot_flag:
        #plot_CAPY_output(body,dataset,output_hyd_namefile[0:-21])
        plot_CAPY_output(body,dataset,'./')

    return True

# convert to WAMIT file from Capytaine generated netCDF file
def convert_CAPYfile_to_WAMITfile(dataset_file,output_hyd_namefile):
    
    dataset = xr.open_dataset(dataset_file,engine='h5netcdf')
    convert_CAPYtoWAMIT_file(dataset,output_hyd_namefile)
    
    return dataset

def plot_CAPY_output(body,dataset,path):
    
    trasl_DoFs = ["Surge","Sway","Heave"]
    rot_DoFs = ["Roll","Pitch","Yaw"]
   
    plt.figure()
    for dof in trasl_DoFs:
        plt.plot(
            dataset['omega']/(np.pi*2),
            np.abs(dataset['diffraction_force'].sel(wave_direction=0, influenced_dof=dof)+
            dataset['Froude_Krylov_force'].sel(wave_direction=0, influenced_dof=dof))/1000,
            label=dof,
            marker='o',
        )
    plt.xlabel('f (Hz)')
    plt.ylabel('Excitation force (kN)')
    plt.legend()
    plt.savefig(path+"\\trasl_exc_force.png") #save as png
    plt.clf()
    
    plt.figure()
    for dof in trasl_DoFs:    
        
        plt.plot(
            dataset['omega']/(np.pi*2),
            (dataset['added_mass'].sel(radiating_dof=dof, influenced_dof=dof))/1000,
            label=dof,
            marker='o',
        )
    plt.xlabel('f (Hz)')
    plt.ylabel('added mass (ton)')
    plt.legend()
    plt.savefig(path+"\\trasl_added_mass.png") #save as png
    plt.clf()
    
    plt.figure()
    for dof in trasl_DoFs:
        plt.plot(
            dataset['omega']/(np.pi*2),
            (dataset['radiation_damping'].sel(radiating_dof=dof, influenced_dof=dof))/1000,
            label=dof,
            marker='o',
        )
    plt.xlabel('f (Hz)')
    plt.ylabel('radiation damping (kN*s/m)')
    plt.legend()
    plt.savefig(path+"\\trasl_radiation_damping.png") #save as png
    plt.clf()
    
    plt.figure()
    for dof in rot_DoFs:
        plt.plot(
            dataset['omega']/(np.pi*2),
            (np.abs(dataset['diffraction_force'].sel(wave_direction=0, influenced_dof=dof)+
            dataset['Froude_Krylov_force'].sel(wave_direction=0, influenced_dof=dof)))/1000,
            label=dof,
            marker='o',
        )
    plt.xlabel('f (Hz)')
    plt.ylabel('Excitation force (kN*m)')
    plt.legend()
    plt.savefig(path+"\\rot_exc_force.png") #save as png
    plt.clf()
    
    plt.figure()
    for dof in rot_DoFs:    
        
        plt.plot(
            dataset['omega']/(np.pi*2),
            (dataset['added_mass'].sel(radiating_dof=dof, influenced_dof=dof))/1000,
            label=dof,
            marker='o',
        )
    plt.xlabel('f (Hz)')
    plt.ylabel('added mass (ton*m^2)')
    plt.legend()
    plt.savefig(path+"\\rot_added_mass.png") #save as png
    plt.clf()
    
    plt.figure()
    for dof in rot_DoFs:
        plt.plot(
            dataset['omega']/(np.pi*2),
            (dataset['radiation_damping'].sel(radiating_dof=dof, influenced_dof=dof))/1000,
            label=dof,
            marker='o',
        )
    plt.xlabel('f (Hz)')
    plt.ylabel('radiation damping (kN*s)')
    plt.legend()
    plt.savefig(path+"\\rot_radiation_damping.png") #save as png
    plt.clf()
    
if __name__ == '__main__':

    # -------- Only this part needs to be changed -----:
    #omega = [0.01, 0.02, 0.04,0.06,0.08,0.1,0.15, 0.2,0.25, 0.3,0.35, 0.4,0.45, 0.5,0.55, 0.6,0.65, 0.7, 0.8,0.9, 1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,4.0,np.infty]
    #omega = [0.01,1,4]
    #omega = [1]
    omega = np.linspace(0.1,4,15)
    #directions = np.linspace(-np.pi,np.pi,8)
    directions = [0]
    dataset = create_hydrodyn_database("prova_Hydraspar_notclipped.msh", omega, directions,show=True)
    isConverted = convert_CAPYtoWAMIT_file(dataset,"prova_Hydraspar")
    #print(dataset)


