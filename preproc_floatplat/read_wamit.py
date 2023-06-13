# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 16:24:00 2022

@author: Giancarlo
"""

import pandas as pd
import numpy as np
import xarray as xr


class WamitData:
    def __init__(self, omega=np.array([]), addedMass=np.array([]), radDamping=np.array([]), excForces=np.array([])):
        self.omega = omega
        self.addedMass = addedMass
        self.radDamping = radDamping
        self.excForces = excForces
        
    def read_wamit(self, filename):
    
        # read radiation coefficients (file .1) --> output a dictionary with omega (rad/s) as key value 
        df1 = pd.read_csv(filename+'.1', skiprows = 0, delim_whitespace=True, na_values=[''], names = ['period', 'i', 'j', 'A', 'B'])    
    
        omega1=np.zeros(len(df1))
        omega1=2*np.pi/df1['period']
        omega1[omega1<0]=0
        omega1=np.unique(omega1)
        omega1=np.sort(omega1)
        self.omega1=omega1
        
        self.addedMass = {}
        self.radDamping = {}
        for kk in range(len(omega1)):
            self.addedMass[omega1[kk]]=np.zeros([6,6])
            self.radDamping[omega1[kk]]=np.zeros([6,6])

        for kk in range(len(df1)):
            if df1['period'].iloc[kk] != -1:
                self.addedMass[2*np.pi/df1['period'].iloc[kk]][df1['i'].iloc[kk]-1,df1['j'].iloc[kk]-1]=df1['A'].iloc[kk]
                self.radDamping[2*np.pi/df1['period'].iloc[kk]][df1['i'].iloc[kk]-1,df1['j'].iloc[kk]-1]=df1['B'].iloc[kk]
            else:
                self.addedMass[0][df1['i'].iloc[kk]-1,df1['j'].iloc[kk]-1]=df1['A'].iloc[kk]
                self.radDamping[0][df1['i'].iloc[kk]-1,df1['j'].iloc[kk]-1]=df1['B'].iloc[kk]
        
        # read excitation force coefficients (file .3) 
        df3 = pd.read_csv(filename+'.3', skiprows = 0, delim_whitespace=True, names = ['period', 'beta', 'i', 'Mod', 'Phase', 'Re', 'Im'])
        omega3=np.zeros(len(df3))
        omega3=2*np.pi/df3['period']
        omega3[omega3<0]=0
        omega3=np.unique(omega3)
        omega3=np.sort(omega3)
        self.omega3=omega3
        
        beta=np.zeros(len(df3))
        beta=df3['beta']
        beta=np.unique(beta)
        beta=np.sort(beta)
        self.beta=beta
        
        self.excForces = {}
        for ii in range(len(omega3)):
            for jj in range(len(beta)):
                self.excForces[(omega3[ii],beta[jj])]=np.zeros([6,2])
            
        for kk in range(len(df3)):
                if df3['period'].iloc[kk] != -1:
                    self.excForces[(2*np.pi/df3['period'].iloc[kk],df3['beta'].iloc[kk])][df3['i'].iloc[kk]-1,0]=df3['Mod'].iloc[kk]
                    self.excForces[(2*np.pi/df3['period'].iloc[kk],df3['beta'].iloc[kk])][df3['i'].iloc[kk]-1,1]=df3['Phase'].iloc[kk]
                else:
                    self.excForces[(0,df3['beta'].iloc[kk])][df3['i'].iloc[kk]-1,0]=df3['Mod'].iloc[kk]
                    self.excForces[(0,df3['beta'].iloc[kk])][df3['i'].iloc[kk]-1,1]=df3['Phase'].iloc[kk]
        
        return True

if __name__=='__main__':
    
    wm=WamitData()
    print()
    
    wm.read_wamit(filename=r'D:\04_Floatech\WP2\float_plat_mod\run_capytaine_SW_rot\SW_Capytaine_rot')
    