import numpy as np

from preproc_floatplat.floatplatgmshes import create_hydraspar_mesh
from preproc_floatplat.floatplatmmhydrost import calc_equil
from preproc_floatplat.func_not_clipped import create_hydraspar_mesh_notclipped, create_hydraspar_brace_mesh_notclipped

from meshmagick import mmio, hydrostatics, mesh

# %% Input data
name_file='test_hydra'


rho_material=7850.
surf_thickness=0.02

# Environment
rho_w=1023
grav=9.81

# Turbine
Turbine_mass_tons=12.5
Turbine_zCoG=17.3

Rated_thrust=19 # kN
hub_height=22 #m

# Central spar / general
draft=9
interface_TwSL=6  # freeboard

spar_central1_radius = 1.0 # m, higher part
spar_central1_height = 7.0
spar_central2_radius = 2.5 # m, lower part
spar_central2_height = 5

heave_plate_radius   = 5.0  # m
heave_plate_thickness= 0.06 # m

rho_ballast=7850

# lateral spar
spar_lateral_N       = 4    # m
spar_lateral_angle   = 50   # deg
spar_lateral_radius  = 1.0  # m
spar_lateral_length  = (draft+interface_TwSL)/np.cos(np.deg2rad(spar_lateral_angle)) # m
alfa_init            = 0    # deg, first spar azimutal angle


# mesh size control
mesh_size_min=0.2 # m
mesh_size_max=0.3 # m

# %%
# create_hydraspar_mesh(name_file, draft, interface_TwSL,
#                          spar_central1_radius, spar_central1_height,
#                          spar_central2_radius, spar_central2_height,
#                          heave_plate_radius, heave_plate_thickness,
#                          spar_lateral_radius, spar_lateral_length,
#                          spar_lateral_N, spar_lateral_angle,
#                          alfa_init,
#                          mesh_size_min, mesh_size_max,
#                          show=True)

# create_hydraspar_mesh_notclipped(name_file, draft, interface_TwSL,
#                          spar_central1_radius, spar_central1_height,
#                          spar_central2_radius, spar_central2_height,
#                          heave_plate_radius, heave_plate_thickness,
#                          spar_lateral_radius, spar_lateral_length,
#                          spar_lateral_N, spar_lateral_angle,
#                          alfa_init,
#                          mesh_size_min, mesh_size_max,
#                          show=True)

create_hydraspar_brace_mesh_notclipped(name_file, draft=9, interface_TwSL=6,
                        spar_central1_radius=1, spar_central1_height=7,
                        spar_central2_radius=2.5, spar_central2_height=5,
                        heave_plate_radius=5, heave_plate_thickness=0.02,
                        spar_lateral_radius=1, spar_lateral_length=23.33585740291,
                        spar_lateral_N=4, spar_lateral_angle=50,
                        brace_height=12, brace_radius=0.5, brace_angle=80,
                        alfa_init=0,
                        mesh_size_min=0.3, mesh_size_max=0.3,
                        show=True)

V, F = mmio.load_MSH(name_file+'.msh')
mymesh = mesh.Mesh(V, F)
# mymesh.show()

Inertia=mymesh.eval_shell_mesh_inertias(rho_medium=rho_material, thickness=surf_thickness)
print('\n\n')
print(f"Structural mass={Inertia.mass:10.1f} kg")
hydraspar_hydrostat = hydrostatics.compute_hydrostatics(mymesh, [0, 0, 0], rho_w, grav)
platform_mass_tons=Inertia.mass/1000
buoyancy_tons=hydraspar_hydrostat['disp_volume']*rho_w/1000
ballast=buoyancy_tons-platform_mass_tons-Turbine_mass_tons
total_mass_tons=platform_mass_tons+Turbine_mass_tons+ballast


h_ballast=ballast/(rho_ballast*np.pi*spar_central2_radius**2)
z_ballast=-draft+heave_plate_thickness+h_ballast/2

COG_plat=Inertia.gravity_center
COG=COG_plat
COG[2]=(COG_plat[2]*platform_mass_tons+Turbine_mass_tons*Turbine_zCoG+ballast*z_ballast)/total_mass_tons
print(f"COG=[{COG[0]:10.3f}, {COG[1]:10.3f}, {COG[2]:10.3f}] m")
hydraspar_hydrostat = hydrostatics.compute_hydrostatics(mymesh, [COG[0], COG[1], COG[2]], rho_w, grav)


print(f"Ballast = {ballast:.2f} tons")

CoG_Z_displ = hydrostatics.displacement_equilibrium(mymesh, total_mass_tons, rho_w, grav, cog=[COG[0], COG[1], COG[2]], reltol=1e-3, verbose=False)

print('\n\n')
print(f'CoG_Z={COG[2]+CoG_Z_displ:.4f} m')
print(f"total displced mass= {hydraspar_hydrostat['disp_mass']:.1f} kg")
K=hydraspar_hydrostat['stiffness_matrix']
print("stiffness matrix ")
print(hydraspar_hydrostat['stiffness_matrix'])

heeling_angle=np.rad2deg(Rated_thrust*1000*(hub_height-COG[2])/K[2,2])
print(f"heeling angle= {heeling_angle:.2f} deg")
print('\n\n')
