# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 14:12:47 2022

@author: Lazzerini Guido
"""

import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# plt.rcParams.update({
    # "text.usetex": False,
    # "font.family": "serif",
    # "font.serif": ["Helvetica"],
# })

plt.rcParams.update({'font.size': 26})
import seaborn as sns
#from plotly.offline import plot
#import plotly.graph_objects as go
#from plotly.subplots import make_subplots
from scipy.stats import pearsonr

#sns.set_theme(style="ticks")

namefile_optimoutput = r"Pop_list_new_SOFTWIND.txt"

#namefile_DEoptimoutput = ".\DE_output.txt"

# Import dataframe with individuals and functions evaluations

df = pd.read_csv(namefile_optimoutput, sep=" ", on_bad_lines='skip')

print(df['fairlead_radius'])
df = df.loc[df['fairlead_radius'] != "ABORT"]

#dfDE = pd.read_csv(namefile_DEoptimoutput, sep=" ", on_bad_lines='skip')

# Define method to calculate obj function and cost functions

# def freq_obj_function_eval(row):
    # obj_TNatYaw = 100
    # return np.abs(obj_TNatYaw-row['TNatYaw'])/obj_TNatYaw
# 
# def obj_function_eval(row):
    # w_freq = 1
    # w_cost = 0
    # obj_TNatYaw = 100
    # return w_freq*np.abs(obj_TNatYaw-row['TNatYaw'])/obj_TNatYaw + w_cost*(np.abs(row['DeltaChainCost'])+np.abs(row['DeltaBracesCost']))

# def chain_length_eval(row):
#     FairleadHeight=8.7
#     WtrDpth = 180
#     fairlead_anchor_distance = np.sqrt((WtrDpth+FairleadHeight)**2+(row['AnchorRadius']-(row['SparDistance']+7.5+21.0))**2)
#     return row['LineLengthFactor']*fairlead_anchor_distance


# def is_in_DE(row):
#     for j in dfDE.index:
#         flagDE = np.abs(row['AnchorRadius']-dfDE['AnchorRadius'][j]) <= 0.01
#     return flagDE

# def color_DE(row):
#     for j in dfDE.index:
#         flagDE = np.abs(row['AnchorRadius']-dfDE['AnchorRadius'][j]) <= 0.1
#         if flagDE:
#             color= 'r'
#         else:
#             color= 'k'
#     return color

# def corrfunc(x, y, ax=None, **kws):
    # """Plot the correlation coefficient in the top left hand corner of a plot."""
    # r, _ = pearsonr(x, y)
    # ax = ax or plt.gca()
    # ax.annotate(r"$\bf{p = " + format(r,'.2f')+"}$", xy=(.1, .9), xycoords=ax.transAxes,weight="bold")
    # 
# Apply methods to create new columns for the dataframe

#df['FreqObjFunc'] = df.apply(lambda row: freq_obj_function_eval(row), axis=1)
#df['ObjFunc'] = df.apply(lambda row: obj_function_eval(row), axis=1)
#df['RoundObjFunc'] = df.apply(lambda row: obj_function_eval(row), axis=1).round(decimals = 2)
#df['ChainLength'] = df.apply(lambda row: chain_length_eval(row), axis=1)
#df['FlagDE'] = df.apply(lambda row: is_in_DE(row), axis=1)
#df['ColorDE'] = df.apply(lambda row: color_DE(row), axis=1)

# counting unique values
#n1 = len(pd.unique(df['MLLF']))
#n2 = len(df['MLLF'])

# print("No.of.unique values of ObjFunc :",n1,"No. of values: ",n2)

# # Simple plot of iteration evolutions

# fig, (ax1, ax2 , ax3,ax4) = plt.subplots(4,1, sharex=True, layout="constrained")
# ax1.plot(dfDE.index+1,dfDE['ObjFunc'],marker='.' ,markersize=10)
# ax1.grid()
# ax1.set_title('(a)\n')
# ax1.set_ylabel(r'ObjFunc (-)')
# ax1ticks = [0.060,0.070,0.080,0.090];
# ax1.set_yticks(ticks = ax1ticks, labels = ['0.060','0.070','0.080','0.090'])

# ax2.plot(dfDE.index+1,dfDE['AnchorRadius'],marker='.' ,markersize=10)
# ax2.grid()
# ax2.set_title('(b)\n')
# ax2.set_ylabel(r'x_{A} (m)')
# ax2ticks = [400,500,600];
# ax2.set_yticks(ticks = ax2ticks, labels = ['400.0','500.0','600.0'])

# ax3.plot(dfDE.index+1,dfDE['LineLengthFactor'],marker='.' ,markersize=10)
# ax3.grid()
# ax3.set_title('(c)\n')
# ax3.set_ylabel(r'MLLF (-)')
# ax3.ylim = [1.07,1.09];
# ax3ticks = [1.075,1.080,1.085];
# ax3.set_yticks(ticks = ax3ticks, labels = ['1.075','1.080','1.085'])

# ax4.plot(dfDE.index+1,dfDE['SparDistance'],marker='.' ,markersize=10)
# ax4.grid()
# ax4.set_xticks(ticks = np.arange(1, 51, step=2))
# ax4ticks = [23.0,24.0,25.0,26.0];
# ax4.set_yticks(ticks = ax4ticks, labels = ['23.00','24.00','25.00','26.00'])

# ax4.set_title('(d)\n')
# ax4.set_xlabel('Iteration (-)')
# ax4.set_ylabel('{x}_{S} (-)')

# fig.get_layout_engine().set(w_pad=2 / 72, h_pad=2 / 72, hspace=0.1,wspace=0.1)
# fig.set_size_inches(20, 20)
# fig.savefig('optimization_181022_0.png', dpi=120)

# Statistical plots with seaborn

#g = sns.pairplot(df, corner=False)
#g.map_lower(sns.kdeplot, levels=3, color=".5")

#df = df[df["YawMax"]<28]
#df = df[df["YawFFTmax"]<7.5]

print(df)
#kind="reg", plot_kws={'line_kws':{'color':'red'}}
t = sns.pairplot(df,diag_kind="hist",\
                  x_vars=["fairlead_radius", "fairlead_height", "MLLF","surge_mean","FFTpeak_Roll","FFTpeak_Pitch","FFTpeak_Yaw"],\
                  y_vars=["fairlead_radius", "fairlead_height", "MLLF","surge_mean","FFTpeak_Roll","FFTpeak_Pitch","FFTpeak_Yaw"])

#t.map(corrfunc)
ax = t.axes[0,0]
ax.yaxis.label.set_color('blue')
ax = t.axes[6,0]
ax.xaxis.label.set_color('blue')
ax.set_visible(True)

ax = t.axes[1,0]
ax.yaxis.label.set_color('blue')
ax = t.axes[6,1]
ax.xaxis.label.set_color('blue')
ax.set_visible(True)

ax = t.axes[2,0]
ax.yaxis.label.set_color('blue')
ax = t.axes[6,2]
ax.xaxis.label.set_color('blue')
ax.set_visible(True)

ax = t.axes[3,0]
ax.yaxis.label.set_color('red')
ax = t.axes[6,3]
ax.xaxis.label.set_color('red')
ax.set_visible(True)

ax = t.axes[4,0]
ax.yaxis.label.set_color('red')
ax = t.axes[6,4]
ax.xaxis.label.set_color('red')
ax.set_visible(True)

ax = t.axes[5,0]
ax.yaxis.label.set_color('red')
ax = t.axes[6,5]
ax.xaxis.label.set_color('red')
ax.set_visible(True)

ax = t.axes[6,0]
ax.yaxis.label.set_color('red')
ax = t.axes[6,6]
ax.xaxis.label.set_color('red')
ax.set_visible(True)

t.fig.set_figheight(40)
t.fig.set_figwidth(40)
t.tight_layout()
t.savefig('DOE_SOFTWIND.png', dpi=900)
plt.show()

#fig, axis = plt.subplots(3,4,figsize=(15,5))

# for j in range(1,4):
#     if(j==1):
#         jv = 'AnchorRadius'
#     elif(j==2):
#         jv = 'LineLengthFactor'
#     elif(j==3):
#         jv = 'SparDistance'
#     #else:
#     #    jv = 'ChainLength' 
#     for i in range(1, 5):
#         if(i == 1):
#                 a = axis[j-1, i-1].scatter(x = 'HeelMean', y = jv,data = df, c = 'YawFFTPeak', cmap = "bone",linewidth=0.1, edgecolor='white')
#                 axis[j-1, i-1].set_xlim([1.5,3.5])
#         elif(i==2):
#                 b = axis[j-1, i-1].scatter(x = 'YawFFTPeak', y = jv,data = df, c = 'YawFFTPeak', cmap = "bone",linewidth=0.1, edgecolor='white')
#                 axis[j-1, i-1].set_xlim([0,15])
#         elif(i==3):
#                 c = axis[j-1, i-1].scatter(x = 'DeltaChainCost', y = jv,data = df, c = 'YawFFTPeak', cmap = "bone",linewidth=0.1, edgecolor='white')
#                 axis[j-1, i-1].set_xlim([-0.30,0.30])        
#         else:
#                 d = axis[j-1, i-1].scatter(x = 'DeltaBracesCost', y = jv,data = df, c = 'YawFFTPeak', cmap = "bone",linewidth=0.1, edgecolor='white')
#                 axis[j-1, i-1].set_xlim([-0.2,0.4])        

# for j in range(1,4):
#     for i in range(1,5):
#         if j == 3:
#             if(i == 1):
#                 axis[j-1, i-1].set_xlabel('HeelMean (deg)')
#                 axis[j-1, i-1].set_xticks(ticks = [1.5,2.0,2.5,3.0,3.5], labels = ['1.5','2.0','2.5','3.0','3.5'])
#             elif(i==2):
#                 axis[j-1, i-1].set_xlabel('YawFFTPeak (deg)')
#                 axis[j-1, i-1].set_xticks(ticks = [0,5,10,15], labels = ['0','5','10','15'])
#             elif(i==3):
#                 axis[j-1, i-1].set_xlabel('DeltaChainCost (-)')
#                 axis[j-1, i-1].set_xticks(ticks = [-0.20,0,0.20], labels = [r'-20\%',r'0\%',r'+20\%'])
#             else:
#                 axis[j-1, i-1].set_xlabel('DeltaBracesCost (-)')
#                 axis[j-1, i-1].set_xticks(ticks = [0,0.2,0.4], labels = ['0\%',r'+20\%',r'+40\%'])
#         else:
#             axis[j-1, i-1].set_xticks([], [])
#         if i == 1:
#             if(j == 1):
#                 axis[j-1, i-1].set_ylabel(r'x_A (m)')
#                 axis[j-1, i-1].set_yticks(ticks = [450,600,750], labels = ['450','600','750'])
#             elif(j==2):
#                 axis[j-1, i-1].set_ylabel(r'MLLF (-)')
#                 axis[j-1, i-1].set_yticks(ticks = [1.04,1.06,1.08], labels = ['1.04','1.06','1.08'])
#             elif(j==3):
#                 axis[j-1, i-1].set_ylabel(r'x_S (m)')
#                 axis[j-1, i-1].set_yticks(ticks = [22,26,30,34], labels = ['22.0','26.0','30.0','34.0'])
#             #else:
#             #     axis[j-1, i-1].set_ylabel(r'l_{moor} (m)')
#             #     axis[j-1, i-1].set_yticks(ticks = [400,500,600,700], labels = ['400.0','500.0','600.0','700.0'])
#         else:
#             axis[j-1,i-1].set_yticks([],[])

# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# fig.colorbar(a, cax=cbar_ax,label='YawFFTPeak')
# fig.set_size_inches(30, 15)
# fig.savefig('optimization_181022_2.png', dpi=200)