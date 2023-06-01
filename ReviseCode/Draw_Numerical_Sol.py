import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Package_MyCode import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=10,fontsize=16)

# DirPath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/Converg_L30_nc3_o3_T_600_8192_HProj_3_Qorder_5_Lap_Dirichlet'
# res = np.load(os.path.join(DirPath,'NumSol_Nspace_360_NT_8192.npy'),allow_pickle=True).item()

DirPath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/Converg_L20_nc2_o3_T_128_4096_HProj_3_Qorder_5_Lap_Dirichlet'
res = np.load(os.path.join(DirPath,'NumSol_Nspace_540_NT_4096.npy'),allow_pickle=True).item()

Z = abs(np.array(res['fval']))
X, Y = np.meshgrid(res['xref'],res['tset'])

# figure = plt.figure()
# axes = figure.add_subplot(111, projection='3d')
# surf = axes.plot_surface(X,Y,Z,cmap = cm.coolwarm)
fig, ax = plt.subplots(figsize=(8,8))

e = ax.imshow(Z, cmap=cm.coolwarm, interpolation='bilinear',
            #   vmin=-4, vmax=4
          extent=[min(res['xref']),max(res['xref']),min(res['tset']),max(res['tset'])])
ax.set_aspect('auto') # you may also use am.imshow(..., aspect="auto") to restore the aspect ratio
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.colorbar(e)
plt.show()

SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
fig.savefig(os.path.join(SavePath,'{}.pdf'.format('NumerSol')),dpi=600,format='pdf')