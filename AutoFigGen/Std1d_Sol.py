import matplotlib.pyplot as plt
import numpy as np
import os
from IO_Func import *

Myplot_Obj = PlotLineStyleObj()
Myplot_Obj.SetColorCycle()
Myplot_Obj.SetLinestyleCycle()
Myplot_Obj.SetPropCycle()
plt.rc('lines', linewidth=3)

CaseName = "../Numerical_Tests/1d_Standard_Soliton/"
DirPath = os.path.join(CaseName,'ME_Conserv/PPC_L40_400_nc3_o3_T_1_20/NumSol_Nspace_400_NT_20.npy')
res = np.load(DirPath,allow_pickle=True).item()
tsets = res['tset']
T_index = [0,int(len(tsets)/2),int(len(tsets)*3/4),-1]

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ALL_Data = np.array(res['fval'])
for index in T_index:
    x = res['xref']
    y = abs(ALL_Data[index,:])
    plt.plot(x,y,label='T='+format(tsets[index],'0.2f'))
ax.grid(True,which="both",ls="--")
ax.legend()
ax.set_xlim([-12,12])
ax.set_xlabel('x')
ax.set_ylabel('|u|')
plt.show()
