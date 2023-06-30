import numpy as np

from preproc_floatplat.floatplatgmshes import create_hydraspar_mesh
from preproc_floatplat.floatplatmmhydrost import inertia_shell_mesh
from meshmagick import hydrostatics, mmio, mesh

# %% Input data
name_file='test_hydra'


rho_material=7850.
surf_thickness=0.03

# Environment
rho_w = 1023 # kg/m^3
grav = 9.81 # m/s^2

# Turbine
Turbine_mass_tons = 1280 # ton
Turbine_zCoG = 97.10 # m

Rated_thrust = 1600 # kN
hub_height = 118.38 #m

# Central spar / general
draft = 40.5 # m
interface_TwSL = 25 # freeboard

spar_central1_radius = 4.5 # m, higher part
spar_central1_height = 31.5
spar_central2_radius = 9 # m, lower part
spar_central2_height = 22.5

heave_plate_radius   = 15.0  # m
heave_plate_thickness= 0.5 # m

rho_ballast=7850

# lateral spar
spar_lateral_N       = 3    # m
spar_lateral_angle   = 53   # deg
spar_lateral_radius  = 4.5  # m
spar_lateral_length  = (draft+interface_TwSL)/np.cos(np.deg2rad(spar_lateral_angle)) # m
alfa_init            = 0    # deg, first spar azimutal angle

brace_height = 30
brace_radius = 2
brace_angle = 80

# mesh size control
mesh_size_min=1 # m
mesh_size_max=1 # m

#  - create_hydraspar_mesh(name_file,
#                         draft, interface_TwSL, 
#                         spar_central1_radius, spar_central1_height,
#                         spar_central2_radius, spar_central2_height,
#                         heave_plate_radius, heave_plate_thickness,
#                         spar_lateral_radius, spar_lateral_length,
#                         spar_lateral_N, spar_lateral_angle,
#                         spar_brace_radius,
#                         spar_brace_angle,spar_brace_height,
#                         alfa_init, 
#                         mesh_size_min, mesh_size_max,
#                         *,show,clip,sealevelboundary):


create_hydraspar_mesh(name_file, 
                        draft, interface_TwSL,
                        spar_central1_radius, spar_central1_height,
                        spar_central2_radius, spar_central2_height,
                        heave_plate_radius, heave_plate_thickness,
                        spar_lateral_radius, spar_lateral_length,
                        spar_lateral_N, spar_lateral_angle,
                        brace_radius,
                        brace_angle,brace_height,
                        alfa_init,
                        mesh_size_min, mesh_size_max,
                        show=True,
                        clip = False)

input_mesh_file = name_file + r".msh"

Inertia = inertia_shell_mesh(input_mesh_file,rho_material,surf_thickness, show = True)

name_file_clipped = 'test_hydra_clipped'

create_hydraspar_mesh(name_file_clipped, 
                        draft, interface_TwSL,
                        spar_central1_radius, spar_central1_height,
                        spar_central2_radius, spar_central2_height,
                        heave_plate_radius, heave_plate_thickness,
                        spar_lateral_radius, spar_lateral_length,
                        spar_lateral_N, spar_lateral_angle,
                        brace_radius,
                        brace_angle,brace_height,
                        alfa_init,
                        mesh_size_min, mesh_size_max,
                        show=True,
                        clip = True)

input_mesh_file_clipped = name_file_clipped + r".msh"

V, F = mmio.load_MSH(input_mesh_file_clipped)
mymesh = mesh.Mesh(V, F)

hydraspar_hydrostat = hydrostatics.compute_hydrostatics(mymesh, [0, 0, 0], rho_w, grav)
platform_mass_tons = Inertia[0].mass/1000
buoyancy_tons = hydraspar_hydrostat['disp_volume']*rho_w/1000
ballast = buoyancy_tons-platform_mass_tons-Turbine_mass_tons
total_mass_tons = platform_mass_tons+Turbine_mass_tons+ballast
z_CoB = hydraspar_hydrostat['buoyancy_center'][2]

h_ballast=ballast/(rho_ballast*np.pi*spar_central2_radius**2)
z_ballast=-draft+heave_plate_thickness+h_ballast/2

print('\n')
print(f"Structural mass={Inertia[0].mass:10.1f} kg")
print(f"Structural mass={(Inertia[0].mass/1000):10.1f} ton")
print('\n')

COG_plat=Inertia[0].gravity_center
COG=COG_plat
COG[2]=(COG_plat[2]*platform_mass_tons+Turbine_mass_tons*Turbine_zCoG+ballast*z_ballast)/total_mass_tons
print(f"COG=[{COG[0]:10.3f}, {COG[1]:10.3f}, {COG[2]:10.3f}] m")
hydraspar_hydrostat = hydrostatics.compute_hydrostatics(mymesh, [COG[0], COG[1], COG[2]], rho_w, grav)
hydraspar_hydrostat_msl = hydrostatics.compute_hydrostatics(mymesh, [0, 0, 0], rho_w, grav)

print(f"Ballast = {ballast:.2f} ton")
print(f"height of ballast = {h_ballast:.2f} m")

CoG_Z_displ = hydrostatics.displacement_equilibrium(mymesh, total_mass_tons, rho_w, grav, cog=[COG[0], COG[1], COG[2]], reltol=1e-3, verbose=False)

z_CoG = COG[2]+CoG_Z_displ

print('\n')
print(f'CoG_Z={COG[2]+CoG_Z_displ:.2f} m')
print(f"total displaced mass= {hydraspar_hydrostat['disp_mass']:.1f} kg")
print(f"total displaced mass= {(hydraspar_hydrostat['disp_mass']/1000):.1f} ton")
print(f"total displaced volume= {hydraspar_hydrostat['disp_volume']:.1f} m^3")



K=hydraspar_hydrostat['stiffness_matrix']
K_msl = hydraspar_hydrostat_msl['stiffness_matrix']

print('\n')
print("stiffness matrix referred to CoG ")
print(K)

print('\n')
print("stiffness matrix referred to [0,0,0] ")
print(K_msl)

print('\n')
heeling_angle=np.rad2deg(Rated_thrust*1000*(hub_height-COG[2])/K[2,2])
print(f"heeling angle= {heeling_angle:.2f} deg (from stiff. matrix referred to CoG)")
print('\n')

C_55 = K_msl[2,2]+hydraspar_hydrostat['disp_volume']*rho_w*grav*z_CoB-hydraspar_hydrostat['disp_mass']*grav*z_CoG
C_buoyancy = + hydraspar_hydrostat['disp_volume']*rho_w*grav*z_CoB
C_mass = -hydraspar_hydrostat['disp_mass']*grav*z_CoG

print(f"C_hydro= {K_msl[2,2]:.1e} N*m/rad ")
print(f"C_buoy= {C_buoyancy:.1e} N*m/rad ")
print(f"C_mass= {C_mass:.1e} N*m/rad ")
print('\n')

print(f"C_55= {C_55:.1f} N*m/rad ")

print('\n')
heeling_angle=np.rad2deg(Rated_thrust*1000*(hub_height)/(C_55))
print(f"heeling angle= {heeling_angle:.2f} deg (from stiff. matrix referred to [0,0,0])")
print('\n')

