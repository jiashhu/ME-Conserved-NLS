import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# DirPath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/Convergence_Res/Spatial_Bi_Soliton_T_Inv/Converg1d_L20_nc3_o1_T_2_100_HProj_2_Qorder_5'
DirPath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/Convergence_Res/Spatial_Bi_Soliton_1d/Converg_L20_nc3_o1_T_2_512_HProj_3_Qorder_5_Lap_Dirichlet/'
res = np.load('/Users/liubocheng/Documents/2022/Code-Collocation_NLS/Convergence_Res/Temporal_Bi_Soliton_T_1d/Converg_L20_nc2_o3_T_2_1024_HProj_1_Qorder_5/NumSol_Nspace_1024_NT_56.npy',allow_pickle=True)
# res = np.load(os.path.join(DirPath,'NumSol_Nspace_1440_NT_512.npy'),allow_pickle=True).item()

Z = abs(np.array(res['fval']))
X, Y = np.meshgrid(res['xref'],res['tset'])

figure = plt.figure()
axes = figure.add_subplot(111, projection='3d')
surf = axes.plot_surface(X,Y,Z,cmap = cm.coolwarm)

plt.show()