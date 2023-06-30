# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 12:29:39 2021

#  Module for processing hydrostatics for FLOATING PLATFORMS of wind turbines with "meshmagick"
#  Run it as an independent script if you need it.
#  Contains: 
#  - functions: write_hydrost_file,
#               get_hydrost,
#               calc_equil,
#               print_mesh_quality,
#               clip_mesh
#
#

@author: Guido Lazzerini

"""
from datetime import datetime


from meshmagick import mmio, hydrostatics, mesh, mesh_clipper, inertia
import os

def write_hydrost_file(input_mesh_file,output_hst_namefile, * ,CoG_Z = 0, rho=1023, g=9.81, show=False):

    V, F = mmio.load_MSH(input_mesh_file)
    mymesh = mesh.Mesh(V, F)
    
    if show:
        temp_show_obj = mesh.Mesh(V, F)
        temp_show_obj.show()
    else:
        None
        
    hydro_cilindro = hydrostatics.compute_hydrostatics(mymesh,[0.0, 0.0, CoG_Z], rho, g)
    
    my_data = [[0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, hydro_cilindro['stiffness_matrix'][0][0]/(g*rho),hydro_cilindro['stiffness_matrix'][0][1]/(g*rho),hydro_cilindro['stiffness_matrix'][0][2]/(g*rho), 0],
            [0, 0,hydro_cilindro['stiffness_matrix'][1][0]/(g*rho),hydro_cilindro['stiffness_matrix'][1][1]/(g*rho),hydro_cilindro['stiffness_matrix'][1][2]/(g*rho), 0],
            [0, 0, hydro_cilindro['stiffness_matrix'][2][0]/(g*rho),hydro_cilindro['stiffness_matrix'][2][1]/(g*rho),hydro_cilindro['stiffness_matrix'][2][2]/(g*rho), 0],
            [0, 0, 0, 0, 0, 0]]
    
    my_data_dim = [[0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, hydro_cilindro['stiffness_matrix'][0][0],hydro_cilindro['stiffness_matrix'][0][1],hydro_cilindro['stiffness_matrix'][0][2], 0],
            [0, 0,hydro_cilindro['stiffness_matrix'][1][0],hydro_cilindro['stiffness_matrix'][1][1],hydro_cilindro['stiffness_matrix'][1][2], 0],
            [0, 0, hydro_cilindro['stiffness_matrix'][2][0],hydro_cilindro['stiffness_matrix'][2][1],hydro_cilindro['stiffness_matrix'][2][2], 0],
            [0, 0, 0, 0, 0, 0]]
    
    K_hst = my_data
    K_dim = my_data_dim

    # Create name of file
    ext_file = ".hst"
    complete_name = "".join([output_hst_namefile,ext_file])
        
    # datetime object containing current date and time
    now = datetime.now()
    
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
    
    complete_name = "".join([output_hst_namefile,dt_string,ext_file]) if os.path.exists(complete_name) else complete_name
    
    file_object = open(complete_name, "w+")

    # skip header
    with file_object:
        for i in range(6):
            for j in range(6):
                line = "   {}   {}   {element:.6E}\n".format(i+1, j+1, element=my_data[i][j])
                file_object.write(line)
    
    return True,complete_name,K_hst, K_dim

def get_hydrost(input_mesh_file, * , CoG_Z= 0, rho=1023, g=9.81,show = False):

    V, F = mmio.load_MSH(input_mesh_file)
    mymesh = mesh.Mesh(V, F)
    
    if show:
        temp_show_obj = mesh.Mesh(V, F)
        temp_show_obj.show()
    else:
        None

    I = mymesh.eval_plain_mesh_inertias(rho)
    
    hydro_cilindro = hydrostatics.compute_hydrostatics(mymesh,[0.0, 0.0, CoG_Z],rho, g)
            
    # Fix cogZ to get the proper I_CM elements when calculating Hydrostatics
    cog_Z = CoG_Z
        
    displaced_mass = hydro_cilindro['disp_mass']
    displaced_vol  = hydro_cilindro['disp_volume']
    
    Ixx = hydro_cilindro['Ixx']
    Iyy = hydro_cilindro['Iyy']
    Izz = hydro_cilindro['Izz']
    
    IxxCM = Ixx - displaced_mass*cog_Z
    IyyCM = Iyy - displaced_mass*cog_Z
    IzzCM = Izz
    
    return displaced_mass, displaced_vol, IxxCM, IyyCM, IzzCM, I

def print_mesh_quality(input_mesh_file, *,show = False,rho = 1025, g = 9.81):
    V, F = mmio.load_MSH(input_mesh_file)
    mymesh = mesh.Mesh(V, F)
    
    if show:
        temp_show_obj = mesh.Mesh(V, F)
        temp_show_obj.show()
    else:
        None
        
    mymesh.print_quality()
    
    return True

def calc_equil(input_mesh_file, cog=[0,0,0], disp_tons=1260,*,rho_water=1023, grav=9.81, reltol=1e-06, verbose=False):
    
    V, F = mmio.load_MSH(input_mesh_file)
    mymesh = mesh.Mesh(V, F)
    
    hydro_cilindro = hydrostatics.compute_hydrostatics(mymesh,cog,rho_water, grav)
    msg_hydro_report=hydrostatics.get_hydrostatic_report(hydro_cilindro)
    print(msg_hydro_report)
    
    # z_corr, rotmat_corr = hydrostatics.full_equilibrium(mymesh, cog, disp_tons, rho_water, grav, reltol=reltol, verbose=verbose)
    # print(f'zcorr={z_corr:.3f} m')
    # print(f'rotmat_corr={rotmat_corr}')
    # return z_corr, rotmat_corr 
    
    z_corr = hydrostatics.displacement_equilibrium(mymesh, disp_tons, rho_water, grav, cog=cog, reltol=reltol, verbose=verbose)
    
    return z_corr

def clip_mesh(input_mesh_file, show=False):
    
    V, F = mmio.load_MSH(input_mesh_file)
    mymesh = mesh.Mesh(V, F)
    
    clipped_mesh = mesh_clipper.MeshClipper(mymesh)

    if show:
        clipped_mesh.lower_mesh.show()
    else:
        None

    return clipped_mesh.lower_mesh
   
def inertia_shell_mesh(input_mesh_file,rho_medium,shell_thickness,*,show=False):
    # Returns list with two inertia objects, one at the mesh origin and the other at the CoG of the body 
    V, F = mmio.load_MSH(input_mesh_file)
    mymesh = mesh.Mesh(V, F)

    if show:
        temp_show_obj = mesh.Mesh(V, F)
        temp_show_obj.show()
    else:
        None

    myBodyInertia_atMeshOrigin = mymesh.eval_shell_mesh_inertias(rho_medium,shell_thickness)
    myBodyInertia_atCoG = myBodyInertia_atMeshOrigin.at_cog

    return myBodyInertia_atMeshOrigin,myBodyInertia_atCoG

if __name__ == '__main__':
    
    isFileWritten, output,K_hst,K_dim = write_hydrost_file('TS_finalconfig2.msh', 'TS_finalconfig2', show = True)
    #m,vol,i1,i2,i3,I = get_hydrost('prova_Hydraspar_notclipped.msh', CoG_Z = 0.0)
    #isQualityPrinted = print_mesh_quality('prova_Hydraspar_notclipped.msh', show = True)
    #clippedMesh = clip_mesh('prova_Hydraspar_notclipped.msh', show = True)
    #calculatedEquil = calc_equil('prova_Hydraspar_notclipped.msh', show = True)
    #inertia = inertia_shell_mesh('TS_finalconfig2.msh',7850,0.01, show = True)

    # Property ".inertia matrix" to have access to the inertia matrix [3x3]
    # Property ".mass" to have access to the mass
    # Property ".gravity_center" to have access to the CoG coordinates [3x1]

    #print(inertia[1].mass)
    #print(inertia[1].gravity_center)
    #print(inertia[1].inertia_matrix)