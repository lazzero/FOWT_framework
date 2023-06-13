# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:23:00 2021

#  float_plat_gmshes.py is Module for creating regular meshes for FLOATING PLATFORMS of wind turbines with GMSH
#  Run it as an independent script if you need it.
#  Contains: 
#  - functions: create_spar_mesh, create_triple_spar_mesh, create_OC3_spar_mesh, create_hydraspar_mesh
#
#

@author: Guido Lazzerini

"""

import gmsh
import sys
import math
import numpy as np

## ------------------------------------------------------------------------------
#  - create_spar_mesh(name_file, spar_radius, spar_height, spar_X, spar_Y, vertical_divisions,
#                  circle_divisions, base_divisions)
## ------------------------------------------------------------------------------

def create_spar_mesh(name_file, spar_radius, spar_height, spar_X, spar_Y, vertical_divisions,circle_divisions, base_divisions,*,show):
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")
    
    
    # We can log all messages for further processing with:
    #gmsh.logger.start()
    
    # This function create a SPAR
    #:::::::::::::z::::::::::::
    #:::::::::::::^::::::::::::
    #:::::::::::::|::::::::::::
    # ----------(-o-)------ ->x
    # ::::::::::(-|-)::::::::::
    # ::::::::::(-|-)::::::::::
    # ::::::::::(-|-)::::::::::
    # ::::::::::(-|-)::::::::::
    # ::::::::::(-|-)::::::::::
    # ::::::::::(-|-)::::::::::
    # :::::::::::_|_:::::::::::
    # ::::::::::::|::::::::::::
    # ::::::::::::|::::::::::::
    
    # Definisco le dimensioni di spar
    
    r_s = spar_radius               #Radius of the spar
    h_s = -spar_height              #Height of the spar
    dx1 = spar_X                 #Distance from center of reference frame to the spar 1 in the X direction
    dy1 = spar_Y                 #Distance from center of reference frame to the spar 1 in the Y direction
    
    # Divisioni archi di cerchio e lati del quadratino centrale
    NN = circle_divisions
    # Divisioni spar in verticale
    NV = vertical_divisions
    # Divisioni base
    NB = base_divisions
    # Fattore di infittimento alle estremità cilindro
    KBV = 0.5
    # Fattore di infittimento alle estremità base
    KB = 0.9
    
    # Creo 5 punti per definire il cerchio superiore della spar
    
    p1 = gmsh.model.geo.addPoint(-r_s+dx1, 0+dy1, 0, tag = 1)
    p2 = gmsh.model.geo.addPoint(0+dx1, r_s+dy1, 0, tag = 2)
    p3 = gmsh.model.geo.addPoint(r_s+dx1, 0+dy1, 0, tag = 3)
    p4 = gmsh.model.geo.addPoint(0+dx1, -r_s+dy1, 0, tag = 4)
    p5 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, 0, tag = 5)
    
    # Creo 5 punti per definire il cerchio inferiore della spar
    
    p6 = gmsh.model.geo.addPoint(-r_s+dx1, 0+dy1, h_s, tag = 6)
    p7 = gmsh.model.geo.addPoint(0+dx1, r_s+dy1, h_s, tag = 7)
    p8 = gmsh.model.geo.addPoint(r_s+dx1, 0+dy1, h_s, tag = 8)
    p9 = gmsh.model.geo.addPoint(0+dx1, -r_s+dy1, h_s, tag = 9)
    p10 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s, tag = 10)
    
    
    # Creo 4 punti per definire il quadratino di base del heave plate
    
    p11 = gmsh.model.geo.addPoint(-r_s*0.4+dx1, 0+dy1, h_s, tag = 11)
    p12 = gmsh.model.geo.addPoint(0+dx1, r_s*0.4+dy1, h_s, tag = 12)
    p13 = gmsh.model.geo.addPoint(r_s*0.4+dx1, 0+dy1, h_s, tag = 13)
    p14 = gmsh.model.geo.addPoint(0+dx1, -r_s*0.4+dy1, h_s, tag = 14)
    
    # Creo archi di cerchio per collegare i punti
    
    ca1 = gmsh.model.geo.addCircleArc(p1,p5,p2)
    ca2 = gmsh.model.geo.addCircleArc(p2,p5,p3)
    ca3 = gmsh.model.geo.addCircleArc(p3,p5,p4)
    ca4 = gmsh.model.geo.addCircleArc(p4,p5,p1)
    
    ca5 = gmsh.model.geo.addCircleArc(p6,p10,p7)
    ca6 = gmsh.model.geo.addCircleArc(p7,p10,p8)
    ca7 = gmsh.model.geo.addCircleArc(p8,p10,p9)
    ca8 = gmsh.model.geo.addCircleArc(p9,p10,p6)
    
    
    # Creo linee verticali
    
    l1 = gmsh.model.geo.addLine(p1,p6)
    l2 = gmsh.model.geo.addLine(p2,p7)
    l3 = gmsh.model.geo.addLine(p3,p8)
    l4 = gmsh.model.geo.addLine(p4,p9)
    
    # Creo linee orizzontali che arrivano al quadratino
    
    l5 = gmsh.model.geo.addLine(p6,p11)
    l6 = gmsh.model.geo.addLine(p7,p12)
    l7 = gmsh.model.geo.addLine(p8,p13)
    l8 = gmsh.model.geo.addLine(p9,p14)
    
    # Creo linee orizzontali del quadratino
    
    l9 = gmsh.model.geo.addLine(p11,p12)
    l10 = gmsh.model.geo.addLine(p12,p13)
    l11 = gmsh.model.geo.addLine(p13,p14)
    l12 = gmsh.model.geo.addLine(p14,p11)
    
    
    # Creo loop per le superfici
    
    cl1 = gmsh.model.geo.addCurveLoop([ca1,l2,-ca5,-l1])
    cl2 = gmsh.model.geo.addCurveLoop([ca2,l3,-ca6,-l2])
    cl3 = gmsh.model.geo.addCurveLoop([ca3,l4,-ca7,-l3])
    cl4 = gmsh.model.geo.addCurveLoop([ca4,l1,-ca8,-l4])
    
    cl5 = gmsh.model.geo.addCurveLoop([ca5,l6,-l9,-l5])
    cl6 = gmsh.model.geo.addCurveLoop([ca6,l7,-l10,-l6])
    cl7 = gmsh.model.geo.addCurveLoop([ca7,l8,-l11,-l7])
    cl8 = gmsh.model.geo.addCurveLoop([ca8,l5,-l12,-l8])
    
    cl9 = gmsh.model.geo.addCurveLoop([l9,l10,l11,l12])
    
    s1 = gmsh.model.geo.addSurfaceFilling([cl1],1)
    s2 = gmsh.model.geo.addSurfaceFilling([cl2],2)
    s3 = gmsh.model.geo.addSurfaceFilling([cl3],3)
    s4 = gmsh.model.geo.addSurfaceFilling([cl4],4)
    
    s5 = gmsh.model.geo.addSurfaceFilling([cl5],5)
    s6 = gmsh.model.geo.addSurfaceFilling([cl6],6)
    s7 = gmsh.model.geo.addSurfaceFilling([cl7],7)
    s8 = gmsh.model.geo.addSurfaceFilling([cl8],8)
    
    s9 = gmsh.model.geo.addSurfaceFilling([cl9],9)
    
    # Sincronizza la Geometria del modello
    gmsh.model.geo.synchronize()
    
    # Per una mesh strutturata mettere True alla variabile "transfinite"
    transfinite = True
    # Creiamo le curve e le superfici transfinite per generare una mesh strutturata
    if transfinite:
        gmsh.model.mesh.setTransfiniteCurve(ca1, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca2, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca3, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca4, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca5, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca6, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca7, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca8, NN)
        gmsh.model.mesh.setTransfiniteCurve(l9, NN)
        gmsh.model.mesh.setTransfiniteCurve(l10, NN)
        gmsh.model.mesh.setTransfiniteCurve(l11, NN)
        gmsh.model.mesh.setTransfiniteCurve(l12, NN)
        gmsh.model.mesh.setTransfiniteCurve(l1, NV,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l2, NV,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l3, NV,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l4, NV,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l5, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l6, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l7, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l8, NB,"Progression",-KB)
    
    
    # Set tutte le superfici transfinite
    for i in range(9):
        gmsh.model.mesh.setTransfiniteSurface(i+1)
    
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    
    # Ricombina i triangoli della mesh in elementi quadrangolari come vuole Nemoh
    gmsh.model.mesh.recombine()
    
    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)
    
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()
        
    gmsh.finalize()
    
    return True

## ------------------------------------------------------------------------------
#  - create_OC3_spar_mesh(name_file, spar_radius1,spar_radius2, spar_height1, spar_height2,spar_draft, spar_X, spar_Y,
#                         vertical_divisions1, vertical_divisions2, vertical_divisions3,
#                         circle_divisions, base_divisions,*,show):
## ------------------------------------------------------------------------------

def create_OC3_spar_mesh(name_file, spar_radius1,spar_radius2, spar_height1, spar_height2,spar_draft, spar_X, spar_Y,
                         vertical_divisions1, vertical_divisions2, vertical_divisions3,
                         circle_divisions, base_divisions,*,show):
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")
    
    
    # We can log all messages for further processing with:
    #gmsh.logger.start()
    
    # This function create a SPAR
    #:::::::::::::z::::::::::::
    #:::::::::::::^::::::::::::
    #:::::::::::::|::::::::::::
    # ----------(-o-)-------->x
    # ::::::::::(-|-):::spar1::
    # ::::::::::(-|-)::::::::::
    # ::::::::::/-|-\::::::::::
    # :::::::::(--|--):::::::::
    # :::::::::(--|--):::::::::
    # :::::::::(--|--):::::::::
    # :::::::::(--|--)::spar2::
    # :::::::::(--|--):::::::::
    # ::::::::::__|__::::::::::
    # :::::::draft|::::::::::::
    # ::::::::::::|::::::::::::
    
    # Definisco le dimensioni di spar
    
    r_s1 = spar_radius1               #Radius of the spar1
    h_s1 = -spar_height1              #Height of the spar1
    r_s2 = spar_radius2
    h_s2 = -(spar_draft-spar_height2)
    h_d = -spar_draft
    dx1 = spar_X                 #Distance from center of reference frame to the spar 1 in the X direction
    dy1 = spar_Y                 #Distance from center of reference frame to the spar 1 in the Y direction
    
    # Divisioni archi di cerchio e lati del quadratino centrale
    NN = circle_divisions
    # Divisioni spar 1 in verticale
    NV1 = vertical_divisions1
    # Divisioni spar 1 - spar 2 in verticale
    NV2 = vertical_divisions2
    # Divisioni spar 3 in verticale
    NV3 = vertical_divisions3
    # Divisioni base
    NB = base_divisions
    # Fattore di infittimento alle estremità cilindro
    KBV = 0.5
    # Fattore di infittimento alle estremità base
    KB = 0.9
    
    # Creo 5 punti per definire il cerchio superiore della spar1
    
    p1 = gmsh.model.geo.addPoint(-r_s1+dx1, 0+dy1, 0, tag = 1)
    p2 = gmsh.model.geo.addPoint(0+dx1, r_s1+dy1, 0, tag = 2)
    p3 = gmsh.model.geo.addPoint(r_s1+dx1, 0+dy1, 0, tag = 3)
    p4 = gmsh.model.geo.addPoint(0+dx1, -r_s1+dy1, 0, tag = 4)
    p5 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, 0, tag = 5)
    
    # Creo 5 punti per definire il cerchio inferiore della spar1
    
    p6 = gmsh.model.geo.addPoint(-r_s1+dx1, 0+dy1, h_s1, tag = 6)
    p7 = gmsh.model.geo.addPoint(0+dx1, r_s1+dy1, h_s1, tag = 7)
    p8 = gmsh.model.geo.addPoint(r_s1+dx1, 0+dy1, h_s1, tag = 8)
    p9 = gmsh.model.geo.addPoint(0+dx1, -r_s1+dy1, h_s1, tag = 9)
    p10 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s1, tag = 10)

    # Creo 5 punti per definire il cerchio superiore della spar2
    
    p11 = gmsh.model.geo.addPoint(-r_s2+dx1, 0+dy1, h_s2, tag = 11)
    p12 = gmsh.model.geo.addPoint(0+dx1, r_s2+dy1, h_s2, tag = 12)
    p13 = gmsh.model.geo.addPoint(r_s2+dx1, 0+dy1, h_s2, tag = 13)
    p14 = gmsh.model.geo.addPoint(0+dx1, -r_s2+dy1, h_s2, tag = 14)
    p15 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s2, tag = 15) 
    
    # Creo 5 punti per definire il cerchio inferiore della spar2
    
    p16 = gmsh.model.geo.addPoint(-r_s2+dx1, 0+dy1, h_d, tag = 16)
    p17 = gmsh.model.geo.addPoint(0+dx1, r_s2+dy1, h_d, tag = 17)
    p18 = gmsh.model.geo.addPoint(r_s2+dx1, 0+dy1, h_d, tag = 18)
    p19 = gmsh.model.geo.addPoint(0+dx1, -r_s2+dy1, h_d, tag = 19)
    p20 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_d, tag = 20) 
    
    # Creo 4 punti per definire il quadratino di base del heave plate
    
    p21 = gmsh.model.geo.addPoint(-r_s2*0.4+dx1, 0+dy1, h_d, tag = 21)
    p22 = gmsh.model.geo.addPoint(0+dx1, r_s2*0.4+dy1, h_d, tag = 22)
    p23 = gmsh.model.geo.addPoint(r_s2*0.4+dx1, 0+dy1, h_d, tag = 23)
    p24 = gmsh.model.geo.addPoint(0+dx1, -r_s2*0.4+dy1, h_d, tag = 24)
    
    # Creo archi di cerchio per collegare i punti
    
    ca1 = gmsh.model.geo.addCircleArc(p1,p5,p2)
    ca2 = gmsh.model.geo.addCircleArc(p2,p5,p3)
    ca3 = gmsh.model.geo.addCircleArc(p3,p5,p4)
    ca4 = gmsh.model.geo.addCircleArc(p4,p5,p1)
    
    ca5 = gmsh.model.geo.addCircleArc(p6,p10,p7)
    ca6 = gmsh.model.geo.addCircleArc(p7,p10,p8)
    ca7 = gmsh.model.geo.addCircleArc(p8,p10,p9)
    ca8 = gmsh.model.geo.addCircleArc(p9,p10,p6)
    
    ca9 = gmsh.model.geo.addCircleArc(p11,p15,p12)
    ca10 = gmsh.model.geo.addCircleArc(p12,p15,p13)
    ca11 = gmsh.model.geo.addCircleArc(p13,p15,p14)
    ca12 = gmsh.model.geo.addCircleArc(p14,p15,p11)
    
    ca13 = gmsh.model.geo.addCircleArc(p16,p20,p17)
    ca14 = gmsh.model.geo.addCircleArc(p17,p20,p18)
    ca15 = gmsh.model.geo.addCircleArc(p18,p20,p19)
    ca16 = gmsh.model.geo.addCircleArc(p19,p20,p16)
    
    # Creo linee verticali spar 1
    
    l1 = gmsh.model.geo.addLine(p1,p6)
    l2 = gmsh.model.geo.addLine(p2,p7)
    l3 = gmsh.model.geo.addLine(p3,p8)
    l4 = gmsh.model.geo.addLine(p4,p9)
    
    # Creo linee verticali tra spar 1 e spar 2
    
    l5 = gmsh.model.geo.addLine(p6,p11)
    l6 = gmsh.model.geo.addLine(p7,p12)
    l7 = gmsh.model.geo.addLine(p8,p13)
    l8 = gmsh.model.geo.addLine(p9,p14)
    
    # Creo linee verticali spar 2
    
    l9 = gmsh.model.geo.addLine(p11,p16)
    l10 = gmsh.model.geo.addLine(p12,p17)
    l11 = gmsh.model.geo.addLine(p13,p18)
    l12 = gmsh.model.geo.addLine(p14,p19)
    
    # Creo linee orizzontali base
    
    l13 = gmsh.model.geo.addLine(p16,p21)
    l14 = gmsh.model.geo.addLine(p17,p22)
    l15 = gmsh.model.geo.addLine(p18,p23)
    l16 = gmsh.model.geo.addLine(p19,p24)
    
    # Creo linee orizzontali del quadratino
    
    l17 = gmsh.model.geo.addLine(p21,p22)
    l18 = gmsh.model.geo.addLine(p22,p23)
    l19 = gmsh.model.geo.addLine(p23,p24)
    l20 = gmsh.model.geo.addLine(p24,p21)
    
    
    # Creo loop per le superfici
    
    cl1 = gmsh.model.geo.addCurveLoop([ca1,l2,-ca5,-l1])
    cl2 = gmsh.model.geo.addCurveLoop([ca2,l3,-ca6,-l2])
    cl3 = gmsh.model.geo.addCurveLoop([ca3,l4,-ca7,-l3])
    cl4 = gmsh.model.geo.addCurveLoop([ca4,l1,-ca8,-l4])
    
    cl5 = gmsh.model.geo.addCurveLoop([ca5,l6,-ca9,-l5])
    cl6 = gmsh.model.geo.addCurveLoop([ca6,l7,-ca10,-l6])
    cl7 = gmsh.model.geo.addCurveLoop([ca7,l8,-ca11,-l7])
    cl8 = gmsh.model.geo.addCurveLoop([ca8,l5,-ca12,-l8])
    
    cl9 = gmsh.model.geo.addCurveLoop([ca9,l10,-ca13,-l9])
    cl10 = gmsh.model.geo.addCurveLoop([ca10,l11,-ca14,-l10])
    cl11 = gmsh.model.geo.addCurveLoop([ca11,l12,-ca15,-l11])
    cl12 = gmsh.model.geo.addCurveLoop([ca12,l9,-ca16,-l12])
    
    cl13 = gmsh.model.geo.addCurveLoop([ca13,l14,-l17,-l13])
    cl14 = gmsh.model.geo.addCurveLoop([ca14,l15,-l18,-l14])
    cl15 = gmsh.model.geo.addCurveLoop([ca15,l16,-l19,-l15])
    cl16 = gmsh.model.geo.addCurveLoop([ca16,l13,-l20,-l16])
    
    cl17 = gmsh.model.geo.addCurveLoop([l17,l18,l19,l20])
    
    s1 = gmsh.model.geo.addSurfaceFilling([cl1],1)
    s2 = gmsh.model.geo.addSurfaceFilling([cl2],2)
    s3 = gmsh.model.geo.addSurfaceFilling([cl3],3)
    s4 = gmsh.model.geo.addSurfaceFilling([cl4],4)
    
    s5 = gmsh.model.geo.addSurfaceFilling([cl5],5)
    s6 = gmsh.model.geo.addSurfaceFilling([cl6],6)
    s7 = gmsh.model.geo.addSurfaceFilling([cl7],7)
    s8 = gmsh.model.geo.addSurfaceFilling([cl8],8)
    
    s9 = gmsh.model.geo.addSurfaceFilling([cl9],9)
    s10 = gmsh.model.geo.addSurfaceFilling([cl10],10)
    s11 = gmsh.model.geo.addSurfaceFilling([cl11],11)
    s12 = gmsh.model.geo.addSurfaceFilling([cl12],12)

    s13 = gmsh.model.geo.addSurfaceFilling([cl13],13)
    s14 = gmsh.model.geo.addSurfaceFilling([cl14],14)
    s15 = gmsh.model.geo.addSurfaceFilling([cl15],15)
    s16 = gmsh.model.geo.addSurfaceFilling([cl16],16)

    s17 = gmsh.model.geo.addSurfaceFilling([cl17],17)


    # Sincronizza la Geometria del modello
    gmsh.model.geo.synchronize()
    
    # Per una mesh strutturata mettere True alla variabile "transfinite"
    transfinite = True
    # Creiamo le curve e le superfici transfinite per generare una mesh strutturata
    if transfinite:
        gmsh.model.mesh.setTransfiniteCurve(ca1, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca2, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca3, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca4, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca5, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca6, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca7, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca8, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca9, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca10, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca11, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca12, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca13, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca14, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca15, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca16, NN)
        gmsh.model.mesh.setTransfiniteCurve(l17, NN)
        gmsh.model.mesh.setTransfiniteCurve(l18, NN)
        gmsh.model.mesh.setTransfiniteCurve(l19, NN)
        gmsh.model.mesh.setTransfiniteCurve(l20, NN)
        gmsh.model.mesh.setTransfiniteCurve(l1, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l2, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l3, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l4, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l5, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l6, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l7, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l8, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l9, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l10, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l11, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l12, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l13, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l14, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l15, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l16, NB,"Progression",-KB)
    
    
    # Set tutte le superfici transfinite
    for i in range(17):
        gmsh.model.mesh.setTransfiniteSurface(i+1)
    
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    
    # Ricombina i triangoli della mesh in elementi quadrangolari come vuole Nemoh
    gmsh.model.mesh.recombine()
    
    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)
        
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()
        
    gmsh.finalize()
    
    return True

## ------------------------------------------------------------------------------
#  - create_triple_spar_mesh(name_file, spar_radius,spar height,hp_radius,hp_thickness,
#                            spar1_X,spar1_Y,spar_vertical_divisions, circle_divisions,
#                            hp_vertical_divisions,hp_sup_divisions,hp_inf_divisions,*,show)
#                  
## ------------------------------------------------------------------------------

def create_triple_spar_mesh(name_file,spar_radius,spar_height,hp_radius,hp_thickness,spar1_X,spar1_Y,spar_vertical_divisions, circle_divisions,hp_vertical_divisions,hp_sup_divisions,hp_inf_divisions,*,show):
    
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")
    
    # Let's build the same model as in `t5.py', but using constructive solid
    # geometry.
    
    # We can log all messages for further processing with:
    #gmsh.logger.start()
    
    # Questo codice crea 3 spar con heave plates nelle posizioni che seguono
    #::::::::::::X:::::::::::
    #::::::::::::|:::::::::::
    #::::::::::::|:::::::::::
    # ::::::::::(1)::::::::::
    #::::::::::/:|:\:::::::::
    #:::::::::/::|::\::::::::
    #Y-------/---O---\-------
    #:::::::/::::|::::\::::::
    #:::::(3)---------(2)::::
    #::::::::::::|:::::::::::
    #::::::::::::|:::::::::::
    #::::::::::::|:::::::::::
    #::::::::::::::::::::::::  
    # Definisco le dimensioni di spar
    
    r_s = spar_radius               #Radius of the spar
    r_hp = hp_radius            #Radius of the heave plate
    h_s =  -spar_height           #Height to the heave plate superior surface
    h_hp = -(spar_height + hp_thickness)          #Height to the heave plate inferior surface
    dx1 =  spar1_X               #Distance from center of reference frame to the spar 1 in the X direction
    dy1 = spar1_Y                 #Distance from center of reference frame to the spar 1 in the Y direction
    dx2 = math.sin(math.radians(30))*(-dx1)  #Distance from center of reference frame to the spar 2 in the X direction (120°)
    dy2 = math.cos(math.radians(30))*(-dx1)  #Distance from center of reference frame to the spar 2 in the Y direction (120°)
    dx3 = dx2               #Distance from center of reference frame to the spar 3 in the X direction
    dy3 = -dy2              #Distance from center of reference frame to the spar 3 in the Y direction
    
    # Creo 5 punti per definire il cerchio superiore della spar
    
    p1 = gmsh.model.geo.addPoint(-r_s+dx1, 0+dy1, 0, tag = 1)
    p2 = gmsh.model.geo.addPoint(0+dx1, r_s+dy1, 0, tag = 2)
    p3 = gmsh.model.geo.addPoint(r_s+dx1, 0+dy1, 0, tag = 3)
    p4 = gmsh.model.geo.addPoint(0+dx1, -r_s+dy1, 0, tag = 4)
    p5 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, 0, tag = 5)
    
    # Creo 5 punti per definire il cerchio inferiore della spar
    
    p6 = gmsh.model.geo.addPoint(-r_s+dx1, 0+dy1, h_s, tag = 6)
    p7 = gmsh.model.geo.addPoint(0+dx1, r_s+dy1, h_s, tag = 7)
    p8 = gmsh.model.geo.addPoint(r_s+dx1, 0+dy1, h_s, tag = 8)
    p9 = gmsh.model.geo.addPoint(0+dx1, -r_s+dy1, h_s, tag = 9)
    p10 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s, tag = 10)
    
    # Creo 5 punti per definire il cerchio superiore del heave plate
    
    p11 = gmsh.model.geo.addPoint(-r_hp+dx1, 0+dy1, h_s, tag = 11)
    p12 = gmsh.model.geo.addPoint(0+dx1, r_hp+dy1, h_s, tag = 12)
    p13 = gmsh.model.geo.addPoint(r_hp+dx1, 0+dy1, h_s, tag = 13)
    p14 = gmsh.model.geo.addPoint(0+dx1, -r_hp+dy1, h_s, tag = 14)
    p15 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s, tag = 15)
    
    # Creo 5 punti per definire il cerchio inferiore del heave plate
    
    p16 = gmsh.model.geo.addPoint(-r_hp+dx1, 0+dy1, h_hp, tag = 16)
    p17 = gmsh.model.geo.addPoint(0+dx1, r_hp+dy1, h_hp, tag = 17)
    p18 = gmsh.model.geo.addPoint(r_hp+dx1, 0+dy1, h_hp, tag = 18)
    p19 = gmsh.model.geo.addPoint(0+dx1, -r_hp+dy1, h_hp, tag = 19)
    p20 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_hp, tag = 20)
    
    # Creo 4 punti per definire il quadratino di base del heave plate
    
    p21 = gmsh.model.geo.addPoint(-r_hp*0.4+dx1, 0+dy1, h_hp, tag = 21)
    p22 = gmsh.model.geo.addPoint(0+dx1, r_hp*0.4+dy1, h_hp, tag = 22)
    p23 = gmsh.model.geo.addPoint(r_hp*0.4+dx1, 0+dy1, h_hp, tag = 23)
    p24 = gmsh.model.geo.addPoint(0+dx1, -r_hp*0.4+dy1, h_hp, tag = 24)
    
    # Creo archi di cerchio per collegare i punti
    
    ca1 = gmsh.model.geo.addCircleArc(p1,p5,p2)
    ca2 = gmsh.model.geo.addCircleArc(p2,p5,p3)
    ca3 = gmsh.model.geo.addCircleArc(p3,p5,p4)
    ca4 = gmsh.model.geo.addCircleArc(p4,p5,p1)
    
    ca5 = gmsh.model.geo.addCircleArc(p6,p10,p7)
    ca6 = gmsh.model.geo.addCircleArc(p7,p10,p8)
    ca7 = gmsh.model.geo.addCircleArc(p8,p10,p9)
    ca8 = gmsh.model.geo.addCircleArc(p9,p10,p6)
    
    ca9 = gmsh.model.geo.addCircleArc(p11,p15,p12)
    ca10 = gmsh.model.geo.addCircleArc(p12,p15,p13)
    ca11 = gmsh.model.geo.addCircleArc(p13,p15,p14)
    ca12 = gmsh.model.geo.addCircleArc(p14,p15,p11)
    
    ca13 = gmsh.model.geo.addCircleArc(p16,p20,p17)
    ca14 = gmsh.model.geo.addCircleArc(p17,p20,p18)
    ca15 = gmsh.model.geo.addCircleArc(p18,p20,p19)
    ca16 = gmsh.model.geo.addCircleArc(p19,p20,p16)
    
    # Creo linee verticali
    
    l1 = gmsh.model.geo.addLine(p1,p6)
    l2 = gmsh.model.geo.addLine(p2,p7)
    l3 = gmsh.model.geo.addLine(p3,p8)
    l4 = gmsh.model.geo.addLine(p4,p9)
    
    # Creo linee orizzontali heave plate superiore
    
    l5 = gmsh.model.geo.addLine(p6,p11)
    l6 = gmsh.model.geo.addLine(p7,p12)
    l7 = gmsh.model.geo.addLine(p8,p13)
    l8 = gmsh.model.geo.addLine(p9,p14)
    
    # Creo linee verticali heave plate
    
    l9 = gmsh.model.geo.addLine(p11,p16)
    l10 = gmsh.model.geo.addLine(p12,p17)
    l11 = gmsh.model.geo.addLine(p13,p18)
    l12 = gmsh.model.geo.addLine(p14,p19)
    
    # Creo linee orizzontali che arrivano al quadratino
    
    l13 = gmsh.model.geo.addLine(p16,p21)
    l14 = gmsh.model.geo.addLine(p17,p22)
    l15 = gmsh.model.geo.addLine(p18,p23)
    l16 = gmsh.model.geo.addLine(p19,p24)
    
    # Creo linee orizzontali del quadratino
    
    l17 = gmsh.model.geo.addLine(p21,p22)
    l18 = gmsh.model.geo.addLine(p22,p23)
    l19 = gmsh.model.geo.addLine(p23,p24)
    l20 = gmsh.model.geo.addLine(p24,p21)
    
    
    # Creo loop per le superfici
    
    cl1 = gmsh.model.geo.addCurveLoop([ca1,l2,-ca5,-l1])
    cl2 = gmsh.model.geo.addCurveLoop([ca2,l3,-ca6,-l2])
    cl3 = gmsh.model.geo.addCurveLoop([ca3,l4,-ca7,-l3])
    cl4 = gmsh.model.geo.addCurveLoop([ca4,l1,-ca8,-l4])
    
    cl5 = gmsh.model.geo.addCurveLoop([ca5,l6,-ca9,-l5])
    cl6 = gmsh.model.geo.addCurveLoop([ca6,l7,-ca10,-l6])
    cl7 = gmsh.model.geo.addCurveLoop([ca7,l8,-ca11,-l7])
    cl8 = gmsh.model.geo.addCurveLoop([ca8,l5,-ca12,-l8])
    
    cl9 = gmsh.model.geo.addCurveLoop([ca9,l10,-ca13,-l9])
    cl10 = gmsh.model.geo.addCurveLoop([ca10,l11,-ca14,-l10])
    cl11 = gmsh.model.geo.addCurveLoop([ca11,l12,-ca15,-l11])
    cl12 = gmsh.model.geo.addCurveLoop([ca12,l9,-ca16,-l12])
    
    cl13 = gmsh.model.geo.addCurveLoop([ca13,l14,-l17,-l13])
    cl14 = gmsh.model.geo.addCurveLoop([ca14,l15,-l18,-l14])
    cl15 = gmsh.model.geo.addCurveLoop([ca15,l16,-l19,-l15])
    cl16 = gmsh.model.geo.addCurveLoop([ca16,l13,-l20,-l16])
    
    cl17 = gmsh.model.geo.addCurveLoop([l17,l18,l19,l20])
    
    s1 = gmsh.model.geo.addSurfaceFilling([cl1],1)
    s2 = gmsh.model.geo.addSurfaceFilling([cl2],2)
    s3 = gmsh.model.geo.addSurfaceFilling([cl3],3)
    s4 = gmsh.model.geo.addSurfaceFilling([cl4],4)
    
    s5 = gmsh.model.geo.addSurfaceFilling([cl5],5)
    s6 = gmsh.model.geo.addSurfaceFilling([cl6],6)
    s7 = gmsh.model.geo.addSurfaceFilling([cl7],7)
    s8 = gmsh.model.geo.addSurfaceFilling([cl8],8)
    
    s9 = gmsh.model.geo.addSurfaceFilling([cl9],9)
    s10 = gmsh.model.geo.addSurfaceFilling([cl10],10)
    s11 = gmsh.model.geo.addSurfaceFilling([cl11],11)
    s12 = gmsh.model.geo.addSurfaceFilling([cl12],12)
    
    s13 = gmsh.model.geo.addSurfaceFilling([cl13],13)
    s14 = gmsh.model.geo.addSurfaceFilling([cl14],14)
    s15 = gmsh.model.geo.addSurfaceFilling([cl15],15)
    s16 = gmsh.model.geo.addSurfaceFilling([cl16],16)
    
    s17 = gmsh.model.geo.addSurfaceFilling([cl17],17)
    
    # Copio e sposto tutte le linee per creare la spar 2 
    
    lines_spar2 = gmsh.model.geo.copy([(1,ca1),(1,ca2),(1,ca3),(1,ca4),(1,ca5),(1,ca6),(1,ca7),(1,ca8),(1,ca9),(1,ca10),(1,ca11),(1,ca12),(1,ca13),(1,ca14),(1,ca15),(1,ca16),(1,l1),(1,l2),(1,l3),(1,l4),(1,l5),(1,l6),(1,l7),(1,l8),(1,l9),(1,l10),(1,l11),(1,l12),(1,l13),(1,l14),(1,l15),(1,l16),(1,l17),(1,l18),(1,l19),(1,l20)])
    
    
    gmsh.model.geo.translate(lines_spar2,-dx1+dx2,-dy1+dy2,0)
    
    
    # Creo i loop dello spar 2
    
    cl1_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[0][1],lines_spar2[17][1],-lines_spar2[4][1],-lines_spar2[16][1]])
    cl2_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[1][1],lines_spar2[18][1],-lines_spar2[5][1],-lines_spar2[17][1]])
    cl3_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[2][1],lines_spar2[19][1],-lines_spar2[6][1],-lines_spar2[18][1]])
    cl4_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[3][1],lines_spar2[16][1],-lines_spar2[7][1],-lines_spar2[19][1]])
    #
    cl5_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[4][1],lines_spar2[21][1],-lines_spar2[8][1],-lines_spar2[20][1]])
    cl6_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[5][1],lines_spar2[22][1],-lines_spar2[9][1],-lines_spar2[21][1]])
    cl7_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[6][1],lines_spar2[23][1],-lines_spar2[10][1],-lines_spar2[22][1]])
    cl8_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[7][1],lines_spar2[20][1],-lines_spar2[11][1],-lines_spar2[23][1]])
    #
    cl9_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[8][1],lines_spar2[25][1],-lines_spar2[12][1],-lines_spar2[24][1]])
    cl10_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[9][1],lines_spar2[26][1],-lines_spar2[13][1],-lines_spar2[25][1]])
    cl11_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[10][1],lines_spar2[27][1],-lines_spar2[14][1],-lines_spar2[26][1]])
    cl12_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[11][1],lines_spar2[24][1],-lines_spar2[15][1],-lines_spar2[27][1]])
    #
    cl13_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[12][1],lines_spar2[29][1],-lines_spar2[32][1],-lines_spar2[28][1]])
    cl14_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[13][1],lines_spar2[30][1],-lines_spar2[33][1],-lines_spar2[29][1]])
    cl15_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[14][1],lines_spar2[31][1],-lines_spar2[34][1],-lines_spar2[30][1]])
    cl16_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[15][1],lines_spar2[28][1],-lines_spar2[35][1],-lines_spar2[31][1]])
    
    cl17_spar2 = gmsh.model.geo.addCurveLoop([lines_spar2[32][1],lines_spar2[33][1],lines_spar2[34][1],lines_spar2[35][1]])
    
    # Creo le superfici dello spar 2
    
    s1_spar2 = gmsh.model.geo.addSurfaceFilling([cl1_spar2],18)
    s2_spar2 = gmsh.model.geo.addSurfaceFilling([cl2_spar2],19)
    s3_spar2 = gmsh.model.geo.addSurfaceFilling([cl3_spar2],20)
    s4_spar2 = gmsh.model.geo.addSurfaceFilling([cl4_spar2],21)
    #
    s5_spar2 = gmsh.model.geo.addSurfaceFilling([cl5_spar2],22)
    s6_spar2 = gmsh.model.geo.addSurfaceFilling([cl6_spar2],23)
    s7_spar2 = gmsh.model.geo.addSurfaceFilling([cl7_spar2],24)
    s8_spar2 = gmsh.model.geo.addSurfaceFilling([cl8_spar2],25)
    #
    s9_spar2 = gmsh.model.geo.addSurfaceFilling([cl9_spar2],26)
    s10_spar2 = gmsh.model.geo.addSurfaceFilling([cl10_spar2],27)
    s11_spar2 = gmsh.model.geo.addSurfaceFilling([cl11_spar2],28)
    s12_spar2 = gmsh.model.geo.addSurfaceFilling([cl12_spar2],29)
    #
    s13_spar2 = gmsh.model.geo.addSurfaceFilling([cl13_spar2],30)
    s14_spar2 = gmsh.model.geo.addSurfaceFilling([cl14_spar2],31)
    s15_spar2 = gmsh.model.geo.addSurfaceFilling([cl15_spar2],32)
    s16_spar2 = gmsh.model.geo.addSurfaceFilling([cl16_spar2],33)
    
    s17_spar2 = gmsh.model.geo.addSurfaceFilling([cl17_spar2],34)
    
    # Copio e sposto tutte le linee per creare la spar 3 
    
    lines_spar3 = gmsh.model.geo.copy([(1,ca1),(1,ca2),(1,ca3),(1,ca4),(1,ca5),(1,ca6),(1,ca7),(1,ca8),(1,ca9),(1,ca10),(1,ca11),(1,ca12),(1,ca13),(1,ca14),(1,ca15),(1,ca16),(1,l1),(1,l2),(1,l3),(1,l4),(1,l5),(1,l6),(1,l7),(1,l8),(1,l9),(1,l10),(1,l11),(1,l12),(1,l13),(1,l14),(1,l15),(1,l16),(1,l17),(1,l18),(1,l19),(1,l20)])
    
    
    gmsh.model.geo.translate(lines_spar3,-dx1+dx3,-dy1+dy3,0)
    
    
    
    # Creo i loop dello spar 3
    
    cl1_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[0][1],lines_spar3[17][1],-lines_spar3[4][1],-lines_spar3[16][1]])
    cl2_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[1][1],lines_spar3[18][1],-lines_spar3[5][1],-lines_spar3[17][1]])
    cl3_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[2][1],lines_spar3[19][1],-lines_spar3[6][1],-lines_spar3[18][1]])
    cl4_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[3][1],lines_spar3[16][1],-lines_spar3[7][1],-lines_spar3[19][1]])
    #
    cl5_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[4][1],lines_spar3[21][1],-lines_spar3[8][1],-lines_spar3[20][1]])
    cl6_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[5][1],lines_spar3[22][1],-lines_spar3[9][1],-lines_spar3[21][1]])
    cl7_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[6][1],lines_spar3[23][1],-lines_spar3[10][1],-lines_spar3[22][1]])
    cl8_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[7][1],lines_spar3[20][1],-lines_spar3[11][1],-lines_spar3[23][1]])
    #
    cl9_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[8][1],lines_spar3[25][1],-lines_spar3[12][1],-lines_spar3[24][1]])
    cl10_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[9][1],lines_spar3[26][1],-lines_spar3[13][1],-lines_spar3[25][1]])
    cl11_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[10][1],lines_spar3[27][1],-lines_spar3[14][1],-lines_spar3[26][1]])
    cl12_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[11][1],lines_spar3[24][1],-lines_spar3[15][1],-lines_spar3[27][1]])
    #
    cl13_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[12][1],lines_spar3[29][1],-lines_spar3[32][1],-lines_spar3[28][1]])
    cl14_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[13][1],lines_spar3[30][1],-lines_spar3[33][1],-lines_spar3[29][1]])
    cl15_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[14][1],lines_spar3[31][1],-lines_spar3[34][1],-lines_spar3[30][1]])
    cl16_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[15][1],lines_spar3[28][1],-lines_spar3[35][1],-lines_spar3[31][1]])
    
    cl17_spar3 = gmsh.model.geo.addCurveLoop([lines_spar3[32][1],lines_spar3[33][1],lines_spar3[34][1],lines_spar3[35][1]])
    
    # Creo le superfici dello spar 3
    
    s1_spar3 = gmsh.model.geo.addSurfaceFilling([cl1_spar3],35)
    s2_spar3 = gmsh.model.geo.addSurfaceFilling([cl2_spar3],36)
    s3_spar3 = gmsh.model.geo.addSurfaceFilling([cl3_spar3],37)
    s4_spar3 = gmsh.model.geo.addSurfaceFilling([cl4_spar3],38)
    #
    s5_spar3 = gmsh.model.geo.addSurfaceFilling([cl5_spar3],39)
    s6_spar3 = gmsh.model.geo.addSurfaceFilling([cl6_spar3],40)
    s7_spar3 = gmsh.model.geo.addSurfaceFilling([cl7_spar3],41)
    s8_spar3 = gmsh.model.geo.addSurfaceFilling([cl8_spar3],42)
    #
    s9_spar3 = gmsh.model.geo.addSurfaceFilling([cl9_spar3],43)
    s10_spar3 = gmsh.model.geo.addSurfaceFilling([cl10_spar3],44)
    s11_spar3 = gmsh.model.geo.addSurfaceFilling([cl11_spar3],45)
    s12_spar3 = gmsh.model.geo.addSurfaceFilling([cl12_spar3],46)
    #
    s13_spar3 = gmsh.model.geo.addSurfaceFilling([cl13_spar3],47)
    s14_spar3 = gmsh.model.geo.addSurfaceFilling([cl14_spar3],48)
    s15_spar3 = gmsh.model.geo.addSurfaceFilling([cl15_spar3],49)
    s16_spar3 = gmsh.model.geo.addSurfaceFilling([cl16_spar3],50)
    
    s17_spar3 = gmsh.model.geo.addSurfaceFilling([cl17_spar3],51)
    
    # Sincronizza la Geometria del modello
    gmsh.model.geo.synchronize()
    
    # Per una mesh strutturata mettere True alla variabile "transfinite"
    transfinite = True
    # Creiamo le curve e le superfici transfinite per generare una mesh strutturata
    if transfinite:
        # Divisioni archi di cerchio e lati del quadratino centrale
        NN = circle_divisions
        # Divisioni spar in verticale
        NV = spar_vertical_divisions
        # Divisioni heave plate in verticale
        NVHV = hp_vertical_divisions
        # Divisioni heave plate superficie superiore
        NHS = hp_sup_divisions
        # Divisioni heave plate superficie inferiore
        NHI = hp_inf_divisions
        gmsh.model.mesh.setTransfiniteCurve(ca1, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca2, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca3, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca4, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca5, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca6, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca7, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca8, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca9, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca10, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca11, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca12, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca13, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca14, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca15, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca16, NN)
        gmsh.model.mesh.setTransfiniteCurve(l17, NN)
        gmsh.model.mesh.setTransfiniteCurve(l18, NN)
        gmsh.model.mesh.setTransfiniteCurve(l19, NN)
        gmsh.model.mesh.setTransfiniteCurve(l20, NN)
        gmsh.model.mesh.setTransfiniteCurve(l1, NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l2, NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l3, NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l4, NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l9, NVHV)
        gmsh.model.mesh.setTransfiniteCurve(l10, NVHV)
        gmsh.model.mesh.setTransfiniteCurve(l11, NVHV)
        gmsh.model.mesh.setTransfiniteCurve(l12, NVHV)
        gmsh.model.mesh.setTransfiniteCurve(l5, NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l6, NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l7, NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l8, NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l13, NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l14, NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l15, NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(l16, NHI,"Bump",0.3)
        for n in range(16):
            gmsh.model.mesh.setTransfiniteCurve(lines_spar2[n][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[16][1], NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[17][1], NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[18][1], NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[19][1], NV,"Bump",0.3)    
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[20][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[21][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[22][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[23][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[24][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[25][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[26][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[27][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[28][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[29][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[30][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[31][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[32][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[33][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[34][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar2[35][1], NN)
        for n in range(16):
            gmsh.model.mesh.setTransfiniteCurve(lines_spar3[n][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[16][1], NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[17][1], NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[18][1], NV,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[19][1], NV,"Bump",0.3)    
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[20][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[21][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[22][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[23][1], NHS,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[24][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[25][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[26][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[27][1], NVHV)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[28][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[29][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[30][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[31][1], NHI,"Bump",0.3)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[32][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[33][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[34][1], NN)
        gmsh.model.mesh.setTransfiniteCurve(lines_spar3[35][1], NN)
    # Set tutte le superfici transfinite
        for i in range(51):
         gmsh.model.mesh.setTransfiniteSurface(i+1)
    
    
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    
    # Ricombina i triangoli della mesh in elementi quadrangolari come vuole Nemoh
    gmsh.model.mesh.recombine()
    
    # ... and save it to disk
    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)
    
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()
    
    
    gmsh.finalize()
    
    return True

## ------------------------------------------------------------------------------
#  - create_hydraspar_mesh(name_file, draft, interface_TwSL, 
#                         spar_central1_radius, spar_central1_height,
#                         spar_central2_radius, spar_central2_height,
#                         heave_plate_radius, heave_plate_thickness,
#                         spar_lateral_radius, spar_lateral_length,
#                         spar_lateral_N, spar_lateral_angle,
#                         alfa_init, 
#                         mesh_size_min, mesh_size_max,
#                         *,show):
## ------------------------------------------------------------------------------

def create_hydraspar_mesh_notclipped(name_file, draft, interface_TwSL, 
                         spar_central1_radius, spar_central1_height,
                         spar_central2_radius, spar_central2_height,
                         heave_plate_radius, heave_plate_thickness,
                         spar_lateral_radius, spar_lateral_length,
                         spar_lateral_N, spar_lateral_angle,
                         alfa_init,
                         mesh_size_min, mesh_size_max,
                         *,show):
    
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")
    
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size_max)
    
    # define first cylinder
    XC1 = 0
    YC1 = 0
    ZC1 = interface_TwSL
    dxC1 = 0
    dyC1 = 0
    dzC1 = -spar_central1_height
    
    # define cone junction
    XCo1 = 0
    YCo1 = 0
    ZCo1 = ZC1 + dzC1
    dxCo1 = 0
    dyCo1 = 0
    dzCo1 = -((draft+interface_TwSL)-(spar_central1_height+spar_central2_height+heave_plate_thickness))
    
    # define second cylinder
    
    XC2 = 0
    YC2 = 0
    ZC2 = ZCo1 + dzCo1
    dxC2 = 0
    dyC2 = 0
    dzC2 = -spar_central2_height

    # define heave plate
    
    XHP1 = 0
    YHP1 = 0
    ZHP1 = ZC2+dzC2
    dxHP1 = 0
    dyHP1 = 0
    dzHP1 = -heave_plate_thickness
    
    RC1 = spar_central1_radius
    RC2 = spar_central2_radius
    RC3 = heave_plate_radius
    
    #define lateral cyilinder
    
    N = spar_lateral_N
    beta = spar_lateral_angle
    l = spar_lateral_length
    alfai = 360/N
    alfa_vec = alfa_init+np.linspace(alfai,360,N)
    
    #gmsh.model.occ.addBox(-R, -R, -R, 2 * R, 2 * R, 2 * R, 1)
    #gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
    C1 = gmsh.model.occ.addCylinder(XC1, YC1, ZC1, dxC1, dyC1, dzC1, RC1, tag = 1)
    Co1 = gmsh.model.occ.addCone(XCo1, YCo1, ZCo1, dxCo1,dyCo1,dzCo1,RC1,RC2, tag = 2)
    C2 = gmsh.model.occ.addCylinder(XC2, YC2, ZC2, dxC2, dyC2, dzC2, RC2, tag = 3)
    HP1 = gmsh.model.occ.addCylinder(XHP1, YHP1, ZHP1, dxHP1, dyHP1, dzHP1, RC3, tag = 4)
    
    Cl = np.empty(N)
    
    for i in range(0,N):
        
        alfa = alfa_vec[i]
        XCl = l*(0.07) *math.sin(math.radians(beta))*math.cos(math.radians(alfa))
        YCl = l*(0.07) *math.sin(math.radians(beta))*math.sin(math.radians(alfa))
        ZCl = -draft + l*(0.07)*math.cos(math.radians(beta))
        RCl = spar_lateral_radius
        dxCl = l*(0.93) *math.sin(math.radians(beta))*math.cos(math.radians(alfa))
        dyCl = l*(0.93) *math.sin(math.radians(beta))*math.sin(math.radians(alfa))
        dzCl = l*(0.93) *math.cos(math.radians(beta))

        Cl = np.append(Cl, [gmsh.model.occ.addCylinder(XCl,YCl,ZCl,dxCl,dyCl,dzCl,RCl, tag = 100+i)], axis=0)
        #gmsh.model.occ.cut([(3,100+i)],[(3,3)],removeTool=False)
        #gmsh.model.occ.remove([(3,100+i)], recursive=False)
        surfaces = gmsh.model.occ.getEntities(2)
        end = len(surfaces)
        print(surfaces)
        if alfa == 180.0:
            gmsh.model.occ.remove([surfaces[(end-3)],surfaces[(end-2)],surfaces[(end-1)]], recursive=True)
        else:
            gmsh.model.occ.remove([surfaces[(end-2)],surfaces[(end-1)]], recursive=True)
        print('last is lateral spar lateral surface --------')
        print(surfaces)
        print('lines --------')
        lines = gmsh.model.occ.getEntities(1)
        print(lines)
        print('---------')

    # Create HP top surface
    gmsh.model.occ.cut([(2, 12)], [(2, 8)], 100)
    
    # Delete volumes and internal surfaces
    gmsh.model.occ.remove([(3,C1),(3,Co1),(3,C2),(3,HP1)], recursive=False)
    gmsh.model.occ.remove([(2,2),(2,3),(2,5),(2,6),(2,8),(2,9),(2,12)], recursive=True)
    print('lines @ this point --------')
    lines = gmsh.model.occ.getEntities(1)
    print(lines)
    print('---------')
    
    # # Create sea level cut plane
    # gmsh.model.occ.addPoint(20, 20, 0, tag = 100)
    # gmsh.model.occ.addPoint(-20, 20, 0, tag = 200)
    # gmsh.model.occ.addPoint(-20, -20, 0, tag = 300)
    # gmsh.model.occ.addPoint(20, -20, 0, tag = 400)
    # gmsh.model.occ.addLine(100, 200, tag = 100)
    # gmsh.model.occ.addLine(300, 200, tag = 200)
    # gmsh.model.occ.addLine(300, 400, tag = 300)
    # gmsh.model.occ.addLine(400, 100, tag = 400)
    # gmsh.model.occ.addCurveLoop([400, 100, -200, 300], 200)
    # gmsh.model.occ.addPlaneSurface([200], 200)
    

    # p = gmsh.model.occ.cut([(2, 1)], [(2, 200)])
    # gmsh.model.occ.remove([p[0][2]],recursive=True)
    
    # for i in range(0,N):
    #     t = gmsh.model.occ.cut([(2, 13+i)], [(2, 200)])
    #     gmsh.model.occ.remove([t[0][2]],recursive=True)
    
    # # remove sea level
    # gmsh.model.occ.remove([(2,200)],recursive=True)
    
    print('lines @ end --------')
    lines = gmsh.model.occ.getEntities(1)
    print(lines)
    print('---------')
    
    # holes = np.empty(N)
    # cut_surf = []
    # k=0
    # for i in range(0,N):
    #     alfa = alfa_vec[i]
    #     if alfa == 180.0:
    #         l1 = gmsh.model.occ.cut([(1, 8)], [(1,404+k+3*i)],removeTool=False)
    #         a = gmsh.model.occ.addCurveLoop([404+k+3*i,l1[0][1][1]])
    #         b = gmsh.model.occ.addCurveLoop([405+k+3*i,l1[0][1][1]])
    #         f1 = gmsh.model.occ.addPlaneSurface([a])
    #         f2 = gmsh.model.occ.addPlaneSurface([b])
    #         print('f1---------')
    #         print(f1)
    #         print('f2---------')
    #         print(f2)
    #         holes[i] = f1
    #         c = gmsh.model.occ.cut([(2, 7)], [(2,f1),(2,f2)],removeObject=False)
    #         print('cut_surf---------')
    #         print(c)
    #         print('---------')
    #         cut_surf.append(c[0])
    #         k=1
    #     else:
    #         a = gmsh.model.occ.addCurveLoop([404+k+3*i])
    #         f = gmsh.model.occ.addPlaneSurface([a])
    #         print('f---------')
    #         print(f)
    #         holes[i] = f
    #         c = gmsh.model.occ.cut([(2, 7)], [(2,f)],removeObject=False)
    #         print('cut_surf---------')
    #         print(c)
    #         print('---------')
    #         cut_surf.append(c[0])


    
    # print('holes---------')
    # print(holes)
    # print('---------')
    
    
    # gmsh.model.occ.remove([(2,7)],recursive=True) 

    # intersect_surf = []
    # intersect_surf.append(cut_surf[0])
    
    # for i in range(0,N-1):
    #     intersect_surf.append(gmsh.model.occ.intersect(intersect_surf[i],cut_surf[i+1])[0])
        
    # create mesh
    gmsh.model.occ.synchronize()    
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()
    gmsh.model.mesh.recombine()

    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)
    
    # show mesh if needed
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()
            
    gmsh.finalize()
        
    return True


def create_OC3_spar_mesh_roty(name_file, spar_radius1,spar_radius2, spar_height1, spar_height2,spar_draft, spar_X, spar_Y,
                         vertical_divisions1, vertical_divisions2, vertical_divisions3,
                         circle_divisions, base_divisions,*,show, rot=None, freeboard=None):
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.add("temporary")
    
    
    # We can log all messages for further processing with:
    #gmsh.logger.start()
    
    # This function create a SPAR
    #:::::::::::::z::::::::::::
    #:::::::::::::^::::::::::::
    #:::::::::::::|::::::::::::
    # ----------(-o-)-------->x
    # ::::::::::(-|-):::spar1::
    # ::::::::::(-|-)::::::::::
    # ::::::::::/-|-\::::::::::
    # :::::::::(--|--):::::::::
    # :::::::::(--|--):::::::::
    # :::::::::(--|--):::::::::
    # :::::::::(--|--)::spar2::
    # :::::::::(--|--):::::::::
    # ::::::::::__|__::::::::::
    # :::::::draft|::::::::::::
    # ::::::::::::|::::::::::::
    
    # Definisco le dimensioni di spar
    
    r_s1 = spar_radius1               #Radius of the spar1
    h_s1 = -spar_height1              #Height of the spar1
    r_s2 = spar_radius2
    h_s2 = -(spar_draft-spar_height2)
    h_d = -spar_draft
    dx1 = spar_X                 #Distance from center of reference frame to the spar 1 in the X direction
    dy1 = spar_Y                 #Distance from center of reference frame to the spar 1 in the Y direction
    
    # Divisioni archi di cerchio e lati del quadratino centrale
    NN = circle_divisions
    # Divisioni spar 1 in verticale
    NV1 = vertical_divisions1
    # Divisioni spar 1 - spar 2 in verticale
    NV2 = vertical_divisions2
    # Divisioni spar 3 in verticale
    NV3 = vertical_divisions3
    # Divisioni base
    NB = base_divisions
    # Fattore di infittimento alle estremità cilindro
    KBV = 0.5
    # Fattore di infittimento alle estremità base
    KB = 0.9
    
    # Creo 5 punti per definire il cerchio superiore della spar1
    if rot!=None:
        rotrad=rot*np.pi/180
    
    if rot==None and freeboard==None:
        # vertcal spar and no freeboard 
        p1 = gmsh.model.geo.addPoint(-r_s1+dx1, 0+dy1, 0, tag = 1)
        p2 = gmsh.model.geo.addPoint(0+dx1, r_s1+dy1, 0, tag = 2)
        p3 = gmsh.model.geo.addPoint(r_s1+dx1, 0+dy1, 0, tag = 3)
        p4 = gmsh.model.geo.addPoint(0+dx1, -r_s1+dy1, 0, tag = 4)
        p5 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, 0, tag = 5)
    elif rot!=None and freeboard==None:
        # rotation with no freeboard
        p1 = gmsh.model.geo.addPoint(-r_s1/np.cos(rotrad)+dx1, 0+dy1, 0, tag = 1)
        p2 = gmsh.model.geo.addPoint(0+dx1, r_s1+dy1, 0, tag = 2)
        p3 = gmsh.model.geo.addPoint(r_s1/np.cos(rotrad)+dx1, 0+dy1, 0, tag = 3)
        p4 = gmsh.model.geo.addPoint(0+dx1, -r_s1+dy1, 0, tag = 4)
        p5 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, 0, tag = 5)
    else:
        # rotation and freeboard (for meshmagick only)
        p1 = gmsh.model.geo.addPoint(-r_s1+dx1, 0+dy1, freeboard, tag = 1)
        p2 = gmsh.model.geo.addPoint(0+dx1, r_s1+dy1, freeboard, tag = 2)
        p3 = gmsh.model.geo.addPoint(r_s1+dx1, 0+dy1, freeboard, tag = 3)
        p4 = gmsh.model.geo.addPoint(0+dx1, -r_s1+dy1, freeboard, tag = 4)
        p5 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, freeboard, tag = 5)
        
    # Creo 5 punti per definire il cerchio inferiore della spar1
    
    p6 = gmsh.model.geo.addPoint(-r_s1+dx1, 0+dy1, h_s1, tag = 6)
    p7 = gmsh.model.geo.addPoint(0+dx1, r_s1+dy1, h_s1, tag = 7)
    p8 = gmsh.model.geo.addPoint(r_s1+dx1, 0+dy1, h_s1, tag = 8)
    p9 = gmsh.model.geo.addPoint(0+dx1, -r_s1+dy1, h_s1, tag = 9)
    p10 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s1, tag = 10)

    # Creo 5 punti per definire il cerchio superiore della spar2
    
    p11 = gmsh.model.geo.addPoint(-r_s2+dx1, 0+dy1, h_s2, tag = 11)
    p12 = gmsh.model.geo.addPoint(0+dx1, r_s2+dy1, h_s2, tag = 12)
    p13 = gmsh.model.geo.addPoint(r_s2+dx1, 0+dy1, h_s2, tag = 13)
    p14 = gmsh.model.geo.addPoint(0+dx1, -r_s2+dy1, h_s2, tag = 14)
    p15 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_s2, tag = 15) 
    
    # Creo 5 punti per definire il cerchio inferiore della spar2
    
    p16 = gmsh.model.geo.addPoint(-r_s2+dx1, 0+dy1, h_d, tag = 16)
    p17 = gmsh.model.geo.addPoint(0+dx1, r_s2+dy1, h_d, tag = 17)
    p18 = gmsh.model.geo.addPoint(r_s2+dx1, 0+dy1, h_d, tag = 18)
    p19 = gmsh.model.geo.addPoint(0+dx1, -r_s2+dy1, h_d, tag = 19)
    p20 = gmsh.model.geo.addPoint(0+dx1, 0+dy1, h_d, tag = 20) 
    
    # Creo 4 punti per definire il quadratino di base del heave plate
    
    p21 = gmsh.model.geo.addPoint(-r_s2*0.4+dx1, 0+dy1, h_d, tag = 21)
    p22 = gmsh.model.geo.addPoint(0+dx1, r_s2*0.4+dy1, h_d, tag = 22)
    p23 = gmsh.model.geo.addPoint(r_s2*0.4+dx1, 0+dy1, h_d, tag = 23)
    p24 = gmsh.model.geo.addPoint(0+dx1, -r_s2*0.4+dy1, h_d, tag = 24)
    
    if rot!=None and freeboard==None:
        gmsh.model.geo.rotate([(0, p6),  (0, p7),  (0, p8),  (0, p9), (0, p10),
                               (0, p11), (0, p12), (0, p13), (0, p14), (0, p15), 
                               (0, p16), (0, p17), (0, p18), (0, p19), (0, p20),
                               (0, p21), (0, p22), (0, p23), (0, p24)], 
                               0, 0, 0, 0, 1, 0, rotrad)
    elif rot!=None and freeboard!=None:
        gmsh.model.geo.rotate([(0, p1),  (0, p2),  (0, p3),  (0, p4), (0, p5),
                               (0, p6),  (0, p7),  (0, p8),  (0, p9), (0, p10),
                               (0, p11), (0, p12), (0, p13), (0, p14), (0, p15), 
                               (0, p16), (0, p17), (0, p18), (0, p19), (0, p20),
                               (0, p21), (0, p22), (0, p23), (0, p24)], 
                               0, 0, 0, 0, 1, 0, rotrad)
    
    
    # Creo archi di cerchio per collegare i punti
    
    if rot==None:
        ca1 = gmsh.model.geo.addCircleArc(p1,p5,p2)
        ca2 = gmsh.model.geo.addCircleArc(p2,p5,p3)
        ca3 = gmsh.model.geo.addCircleArc(p3,p5,p4)
        ca4 = gmsh.model.geo.addCircleArc(p4,p5,p1)
    else:
        ca1 = gmsh.model.geo.addEllipseArc(p1,p5,p2,p2)
        ca2 = gmsh.model.geo.addEllipseArc(p2,p5,p3,p3)
        ca3 = gmsh.model.geo.addEllipseArc(p3,p5,p4,p4)
        ca4 = gmsh.model.geo.addEllipseArc(p4,p5,p1,p1)
        
    ca5 = gmsh.model.geo.addCircleArc(p6,p10,p7)
    ca6 = gmsh.model.geo.addCircleArc(p7,p10,p8)
    ca7 = gmsh.model.geo.addCircleArc(p8,p10,p9)
    ca8 = gmsh.model.geo.addCircleArc(p9,p10,p6)
    
    ca9 = gmsh.model.geo.addCircleArc(p11,p15,p12)
    ca10 = gmsh.model.geo.addCircleArc(p12,p15,p13)
    ca11 = gmsh.model.geo.addCircleArc(p13,p15,p14)
    ca12 = gmsh.model.geo.addCircleArc(p14,p15,p11)
    
    ca13 = gmsh.model.geo.addCircleArc(p16,p20,p17)
    ca14 = gmsh.model.geo.addCircleArc(p17,p20,p18)
    ca15 = gmsh.model.geo.addCircleArc(p18,p20,p19)
    ca16 = gmsh.model.geo.addCircleArc(p19,p20,p16)
    
    # Creo linee verticali spar 1
    
    l1 = gmsh.model.geo.addLine(p1,p6)
    l2 = gmsh.model.geo.addLine(p2,p7)
    l3 = gmsh.model.geo.addLine(p3,p8)
    l4 = gmsh.model.geo.addLine(p4,p9)
    
    # Creo linee verticali tra spar 1 e spar 2
    
    l5 = gmsh.model.geo.addLine(p6,p11)
    l6 = gmsh.model.geo.addLine(p7,p12)
    l7 = gmsh.model.geo.addLine(p8,p13)
    l8 = gmsh.model.geo.addLine(p9,p14)
    
    # Creo linee verticali spar 2
    
    l9 = gmsh.model.geo.addLine(p11,p16)
    l10 = gmsh.model.geo.addLine(p12,p17)
    l11 = gmsh.model.geo.addLine(p13,p18)
    l12 = gmsh.model.geo.addLine(p14,p19)
    
    # Creo linee orizzontali base
    
    l13 = gmsh.model.geo.addLine(p16,p21)
    l14 = gmsh.model.geo.addLine(p17,p22)
    l15 = gmsh.model.geo.addLine(p18,p23)
    l16 = gmsh.model.geo.addLine(p19,p24)
    
    # Creo linee orizzontali del quadratino
    
    l17 = gmsh.model.geo.addLine(p21,p22)
    l18 = gmsh.model.geo.addLine(p22,p23)
    l19 = gmsh.model.geo.addLine(p23,p24)
    l20 = gmsh.model.geo.addLine(p24,p21)
    
    
    # Creo loop per le superfici
    
    cl1 = gmsh.model.geo.addCurveLoop([ca1,l2,-ca5,-l1])
    cl2 = gmsh.model.geo.addCurveLoop([ca2,l3,-ca6,-l2])
    cl3 = gmsh.model.geo.addCurveLoop([ca3,l4,-ca7,-l3])
    cl4 = gmsh.model.geo.addCurveLoop([ca4,l1,-ca8,-l4])
    
    cl5 = gmsh.model.geo.addCurveLoop([ca5,l6,-ca9,-l5])
    cl6 = gmsh.model.geo.addCurveLoop([ca6,l7,-ca10,-l6])
    cl7 = gmsh.model.geo.addCurveLoop([ca7,l8,-ca11,-l7])
    cl8 = gmsh.model.geo.addCurveLoop([ca8,l5,-ca12,-l8])
    
    cl9 = gmsh.model.geo.addCurveLoop([ca9,l10,-ca13,-l9])
    cl10 = gmsh.model.geo.addCurveLoop([ca10,l11,-ca14,-l10])
    cl11 = gmsh.model.geo.addCurveLoop([ca11,l12,-ca15,-l11])
    cl12 = gmsh.model.geo.addCurveLoop([ca12,l9,-ca16,-l12])
    
    cl13 = gmsh.model.geo.addCurveLoop([ca13,l14,-l17,-l13])
    cl14 = gmsh.model.geo.addCurveLoop([ca14,l15,-l18,-l14])
    cl15 = gmsh.model.geo.addCurveLoop([ca15,l16,-l19,-l15])
    cl16 = gmsh.model.geo.addCurveLoop([ca16,l13,-l20,-l16])
    
    cl17 = gmsh.model.geo.addCurveLoop([l17,l18,l19,l20])

    
    s1 = gmsh.model.geo.addSurfaceFilling([cl1],1)
    s2 = gmsh.model.geo.addSurfaceFilling([cl2],2)
    s3 = gmsh.model.geo.addSurfaceFilling([cl3],3)
    s4 = gmsh.model.geo.addSurfaceFilling([cl4],4)
    
    s5 = gmsh.model.geo.addSurfaceFilling([cl5],5)
    s6 = gmsh.model.geo.addSurfaceFilling([cl6],6)
    s7 = gmsh.model.geo.addSurfaceFilling([cl7],7)
    s8 = gmsh.model.geo.addSurfaceFilling([cl8],8)
    
    s9 = gmsh.model.geo.addSurfaceFilling([cl9],9)
    s10 = gmsh.model.geo.addSurfaceFilling([cl10],10)
    s11 = gmsh.model.geo.addSurfaceFilling([cl11],11)
    s12 = gmsh.model.geo.addSurfaceFilling([cl12],12)

    s13 = gmsh.model.geo.addSurfaceFilling([cl13],13)
    s14 = gmsh.model.geo.addSurfaceFilling([cl14],14)
    s15 = gmsh.model.geo.addSurfaceFilling([cl15],15)
    s16 = gmsh.model.geo.addSurfaceFilling([cl16],16)

    s17 = gmsh.model.geo.addSurfaceFilling([cl17],17)


                
    # Sincronizza la Geometria del modello
    gmsh.model.geo.synchronize()

    
    # Per una mesh strutturata mettere True alla variabile "transfinite"
    transfinite = True
    # Creiamo le curve e le superfici transfinite per generare una mesh strutturata
    if transfinite:
        gmsh.model.mesh.setTransfiniteCurve(ca1, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca2, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca3, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca4, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca5, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca6, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca7, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca8, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca9, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca10, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca11, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca12, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca13, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca14, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca15, NN)
        gmsh.model.mesh.setTransfiniteCurve(ca16, NN)
        gmsh.model.mesh.setTransfiniteCurve(l17, NN)
        gmsh.model.mesh.setTransfiniteCurve(l18, NN)
        gmsh.model.mesh.setTransfiniteCurve(l19, NN)
        gmsh.model.mesh.setTransfiniteCurve(l20, NN)
        gmsh.model.mesh.setTransfiniteCurve(l1, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l2, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l3, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l4, NV1,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l5, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l6, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l7, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l8, NV2,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l9, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l10, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l11, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l12, NV3,"Bump",KBV)
        gmsh.model.mesh.setTransfiniteCurve(l13, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l14, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l15, NB,"Progression",-KB)
        gmsh.model.mesh.setTransfiniteCurve(l16, NB,"Progression",-KB)
    
    
        # Set tutte le superfici transfinite
        for i in range(17):
            gmsh.model.mesh.setTransfiniteSurface(i+1)
    
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    
    # Ricombina i triangoli della mesh in elementi quadrangolari come vuole Nemoh
    gmsh.model.mesh.recombine()
    
    # Create name of file
    ext_file = ".msh"
    complete_name = "".join([name_file,ext_file])
    # ... and save it to disk
    gmsh.write(complete_name)
        
    if '-nopopup' not in sys.argv:
        if show:
            gmsh.fltk.run()
        
    gmsh.finalize()
    
    return True


if __name__ == '__main__':
    
    #  - create_triple_spar_mesh(name_file, spar_radius, spar_height,hp_radius,hp_thickness,
    #                            spar1_X,spar1_Y,spar_vertical_divisions, circle_divisions,
    #                            hp_vertical_divisions,hp_sup_divisions,hp_inf_divisions,*,show)
    
    #isMeshCreated = create_spar_mesh('prova_creazione', 7.5, 120, 0, 0, 80,8, 10, show = True)
    #isTSMeshCreated = create_triple_spar_mesh("TS_finalconfig2", 7.5, 53.964, 11.25, 0.5,26, 0, 40, 10, 3, 8, 7, show = True)
    #isOC3MeshCreated = create_OC3_spar_mesh('OC3_2209', 6.5/2, 9.4/2, 4, 108, 120, 0, 0, 10, 15, 80, 8, 10,show= True)

    isHydraSparMeshCreated = create_hydraspar_mesh_notclipped('prova_Hydraspar_notclipped', 9.0, 6.0, 
                              1.0, 7.0,
                              2.5, 5.0,
                              5.0, 0.06,
                              1, 30.0,
                              3, 50.0,
                              0.0,
                              0.3, 0.3,
                              show = True)