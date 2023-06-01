import os
from Package_MyCode import FO
import numpy as np
import matplotlib.pyplot as plt
from Package_MyCode import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=10,fontsize=16)

FoldName1 = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Temporal_Standard_Soliton_1d/Converg_L20_nc2_o3_T_1_1024_HProj_1_Qorder_5'
FoldName2 = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Temporal_Standard_Soliton_1d/NoP_Converg_L20_nc2_o3_T_1_1024_HProj_1_Qorder_5'
def GetEEvsDt(FoldName): 
    t_set = []
    ee_set = []
    for err_file_name in FO.Get_File_List(FoldName):
        if err_file_name.startswith('dt_'):
            data = np.load(os.path.join(FoldName,err_file_name),allow_pickle=True).item()
            t_set.append(data['tau'])
            energyerr = np.array([singleengy-data['endt_energ'][0] for singleengy in data['endt_energ']])
            ee_set.append(max(abs(energyerr)))
    indice = np.argsort(np.array(t_set)).copy()
    t_set = np.sort(np.array(t_set))
    newarr = np.array(ee_set)
    newarr = newarr[indice]
    return t_set, newarr

t_set1, ee_set1 = GetEEvsDt(FoldName1)
t_set2, ee_set2 = GetEEvsDt(FoldName2)
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.loglog(t_set1,ee_set1,'o-',label='The proposed method')
ax.loglog(t_set2,ee_set2,'^--',label='The standard collcation method')
plt.grid(True,ls='--')
ax.set_xlabel('$\\tau$')
ax.set_ylabel('Energy error')
ax.legend(loc='upper left')
plt.show()

SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
fig.savefig(os.path.join(SavePath,'{}.pdf'.format('Example1EEvsdt')),dpi=600,format='pdf')
