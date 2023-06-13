# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 19:39:46 2022

"""

import sys
sys.path.append(r"D:\GitHub\FOWT_optim_test\preproc_floatplat")

from floatplatgmshes import create_OC3_spar_mesh
from floatplatmmhydrost import *
from floatplatcapyhydrodyn import *

import os
import numpy as np

SW_gmsh_file=r'./run_capytaine_SW/SW_gmsh_file'
SW_hydro=r'./run_capytaine_SW/SW_hydro'

resolution_idx=0.5

# generate geometry
create_OC3_spar_mesh(SW_gmsh_file, spar_radius1=5.6, spar_radius2=9, 
                                   spar_height1=5.4, spar_height2=78,
                                   spar_draft=91.4, 
                                   spar_X=0, spar_Y=0,
                                   vertical_divisions1=int(12*resolution_idx), 
                                   vertical_divisions2=int(12*resolution_idx), 
                                   vertical_divisions3=int(50*resolution_idx),
                                   circle_divisions=int(10*resolution_idx), 
                                   base_divisions=int(10*resolution_idx),
                                   show=True)
# run hydrostatic calculations
if os.path.exists(SW_hydro+".hst"):
    os.remove(SW_hydro+".hst")
    
write_hydrost_file(SW_gmsh_file+".msh",SW_hydro, CoG_Z = 0, rho=1023, g=9.81, show=True)

# run hydrodynamic calculations
omegas=np.linspace(0.05,5,5)
omegas=np.append(omegas,10.0)

if os.path.exists(SW_hydro+".1"):
    os.remove(SW_hydro+".1")
if os.path.exists(SW_hydro+".3"):
    os.remove(SW_hydro+".3")

wave_directions=np.linspace(-np.pi,np.pi,3)
data = create_hydrodyn_database(SW_gmsh_file+".msh", omegas, wave_directions, show=False, depth = 200.0)
convert_CAPYtoWAMIT_file(data,SW_hydro)


