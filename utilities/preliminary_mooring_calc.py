# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 22:40:22 2022

@author: Giancarlo
"""
### mooring preliminary design ###

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# water depth
depth=100.0 #m
rho_w=1025.0
g=9.81

# Safety factor
FS=1.0

FH_max=1650.0 #kN

# R3 stud link chain
dc=0.0 # chain diameter
w=(0.0219*dc**2*9.81/1000-(rho_w*g*np.pi*(dc/1000)**2/4)/1000)
Q=0.0
Q_vec=np.array([])
dc_vec=np.array([])
while Q<FH_max+w*depth:
    dc+=1.0
    Q_w=103.84*(44 - 0.08*dc)
    w=(0.0219*dc**2*9.81/1000-(rho_w*g*np.pi*(dc/1000)**2/4)/1000)
    Q=Q_w*w
    Q_vec=np.append(Q_vec,Q)
    dc_vec=np.append(dc_vec,dc)

print(f'\n\n diam={dc:.1f} mm  w={w*1000/9.81:.2f} kg/m Tension allowable = {Q:.2f} kN')

#Ee=120.28 - 0.53*dc

# line angle with respect to elongation direction
alpha=0.0 #deg

# admissible elongation
Dx_adm=25.0

# scope of the mooring S/d = line length/depth
S_d=np.sqrt(2*Q/(w*depth)-1) # 

# Mooring line length
S=S_d*depth
print(f'\n\n Minimum mooring line length: S={S:.2f} m')
# Mooring line footprint at FH_max
a=(Q-w*depth)/w
X_max=a*np.arcsinh(S/a)
print(f' Maximum mooring line footprint at FH_max: Xmax={X_max:.2f} m')

# Mooring line footprint and horizontal force at the mooring system rest position
X_0=X_max-Dx_adm
print(f'\n Mooring line footprint at the mooring system rest position: Xmax={X_0:.2f} m')
X_func = lambda a : S-depth*np.sqrt(1+2*a/depth)+a*np.arccosh(1+depth/a)
X_X0_func=lambda a : X_func(a)-X_0
plt.plot(X_func(np.linspace(1.0,5000,500)),w*np.linspace(1.0,5000,500))
plt.axhline(y=Q+w*depth,color='r')
plt.plot(np.array([X_max,X_max]),np.array([0.0, Q+w*depth]),'r')
plt.grid('on')
plt.xlabel('X (m)')
plt.ylabel('T_H (N)')

a_0=brentq(X_X0_func,0.1,10000)
TH_0=a_0*w
print(f' Mooring line horizontal force at the mooring system rest position: TH_0={TH_0:.2f} kN')
plt.axhline(y=TH_0,color='g')
plt.plot(np.array([X_0,X_0]),np.array([0.0, TH_0]),'g')



