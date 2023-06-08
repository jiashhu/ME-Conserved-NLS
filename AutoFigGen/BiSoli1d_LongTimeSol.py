import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import cm
from IO_Func import PlotLineStyleObj

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=10,fontsize=16)

DirPath = '../Numerical_Tests/1d_Bi_Soliton/LongTSol/PPC_L20_640_nc2_o3_T_128_4096_0'
res = np.load(os.path.join(DirPath,'NumSol_Nspace_640_NT_4096.npy'),allow_pickle=True).item()

Z = abs(np.array(res['fval']))
X, Y = np.meshgrid(res['xref'],res['tset'])
fig, ax = plt.subplots(figsize=(8,8))

e = ax.imshow(Z, cmap=cm.coolwarm, interpolation='bilinear',
            #   vmin=-4, vmax=4
          extent=[min(res['xref']),max(res['xref']),min(res['tset']),max(res['tset'])])
ax.set_aspect('auto') # you may also use am.imshow(..., aspect="auto") to restore the aspect ratio
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.colorbar(e)
plt.show()

# SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
# fig.savefig(os.path.join(SavePath,'{}.pdf'.format('NumerSol')),dpi=600,format='pdf')