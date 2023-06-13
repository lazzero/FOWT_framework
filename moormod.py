# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 22:57:18 2021

@author: Giancarlo Troise & Guido Lazzerini
"""
import numpy as np

from pyFAST.input_output.fast_input_file import FASTInputFile
from pyQBlade.qblade_input_file import QBladeInputFile

def moor_config_openfast(xx,
                         turbModel,
                         filepath='unused', filepath_mod='unused'):

        # Check if "Moordyn" file was provided
        fst_data = FASTInputFile(filepath_mod+"\\"+turbModel['FSTMODFILENAME'])
        if fst_data['CompMooring'] == 0:
            return -1

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
        
        # read hydrodyn input template file
        HD_Data = FASTInputFile(filepath+'\\'+turbModel['HYDFILENAME'])
        WtrDpth = HD_Data['WtrDpth']
        AnchorDepth = WtrDpth

        # read moordyn input template file
        MD_data = FASTInputFile(filepath+'\\'+turbModel['MRDFILENAME'])
        MD_data.NLines=LineNumber
        
        for kk in range(LineNumber):
            # set fairlead and anchor postions
            MD_data['ConnectionProp'][kk,2]=AnchorRadius*np.cos(np.deg2rad(kk*360/LineNumber))
            MD_data['ConnectionProp'][kk,3]=AnchorRadius*np.sin(np.deg2rad(kk*360/LineNumber))
            MD_data['ConnectionProp'][kk,4]=-AnchorDepth
            
            MD_data['ConnectionProp'][kk+LineNumber,2]=FairleadRadius*np.cos(np.deg2rad(kk*360/LineNumber))
            MD_data['ConnectionProp'][kk+LineNumber,3]=FairleadRadius*np.sin(np.deg2rad(kk*360/LineNumber))
            MD_data['ConnectionProp'][kk+LineNumber,4]=FairleadHeight
        
            # set mooring line length
            MD_data['LineProp'][kk,2]=LineLengthFactor*np.sqrt((MD_data['ConnectionProp'][kk,2]-MD_data['ConnectionProp'][kk+LineNumber,2])**2+ \
                                                               (MD_data['ConnectionProp'][kk,3]-MD_data['ConnectionProp'][kk+LineNumber,3])**2+ \
                                                               (MD_data['ConnectionProp'][kk,4]-MD_data['ConnectionProp'][kk+LineNumber,4])**2)
            
            
        # write mooring file
        MD_data.write(filepath_mod+'\\'+turbModel['MRDMODFILENAME'])
        out_data = MD_data
                
        return out_data

def moor_config_qblade(xx,
                       turbModel,
                       filepath='unused',filepath_mod='unused'):
    
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

    sub_Data = QBladeInputFile(filepath +"\\"+ turbModel['SUBFILENAME'])
    WtrDpth=sub_Data['WATERDEPTH']
    AnchorDepth = WtrDpth
    
    for kk in range(LineNumber):
        XA = AnchorRadius*np.cos(np.deg2rad(kk*360/LineNumber))
        YA = AnchorRadius*np.sin(np.deg2rad(kk*360/LineNumber))
        ZA = -AnchorDepth
        XF = FairleadRadius*np.cos(np.deg2rad(kk*360/LineNumber))
        YF = FairleadRadius*np.sin(np.deg2rad(kk*360/LineNumber))
        ZF = FairleadHeight
        
        anchor_pos= "GRD_%.2f_%.2f" % (XA, YA)
        floating_pos = "FLT_%.2f_%.2f_%.2f" % (XF, YF, ZF)
        
        sub_Data.modTabElm('MOORMEMBERS',kk,'CONN_1',floating_pos)
        sub_Data.modTabElm('MOORMEMBERS',kk,'CONN_2',anchor_pos)
        sub_Data.modTabElm('MOORMEMBERS',kk,'Len.[m]',LineLengthFactor*np.sqrt((XA-XF)**2+(YA-YF)**2+(ZA-ZF)**2))
        
    sub_Data.write(filepath_mod +"\\"+ turbModel['SUBMODFILENAME'])
    out_data = sub_Data
    
    return out_data