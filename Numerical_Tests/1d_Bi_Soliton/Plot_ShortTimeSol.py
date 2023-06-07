import matplotlib.pyplot as plt
import numpy as np
import os
from Package_MyCode import PlotLineStyleObj, FO

Myplot_Obj = PlotLineStyleObj()
Myplot_Obj.SetColorCycle()
Myplot_Obj.SetLinestyleCycle()
Myplot_Obj.SetPropCycle()
plt.rc('lines', linewidth=3)

DirPath = './Temp_Conv/Converg_L20_nc4_o3_T_2_2048_HProj_1/NumSol_Nspace_2048_NT_20.npy'
res = np.load(DirPath,allow_pickle=True).item()
tsets = res['tset']
T_index = [0,int(len(tsets)/3),int(len(tsets)*2/3),-1]

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ALL_Data = np.array(res['fval'])
for index in T_index:
    x = res['xref']
    y = abs(ALL_Data[index,:])
    plt.plot(x,y,label='T='+format(tsets[index],'0.2f'))
ax.grid(True,which="both",ls="--")
ax.legend()
ax.set_xlim([-8,8])
ax.set_xlabel('x')
ax.set_ylabel('|u|')
plt.show()

# SavePath = '/Users/hjs/Documents/2022/Collocation_NLS/Draw_Pic/Pic'
# fig.savefig(os.path.join(SavePath,'{}.eps'.format('NumSim_Bi_Soliton_T_1d')),dpi=600,format='eps')