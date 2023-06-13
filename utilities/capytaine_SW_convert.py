# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:09:47 2022

@author: Giancarlo
"""

from floatplatgmshes import create_OC3_spar_mesh
from floatplatmmhydrost import *
from floatplatcapyhydrodyn import *

dataset_file='./run_capytaine_SW/SW_gmsh_file.cpt'
output_hyd_namefile='./run_capytaine_SW/SW_hydro'
dataset = convert_CAPYfile_to_WAMITfile(dataset_file,output_hyd_namefile)
