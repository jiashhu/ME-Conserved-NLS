import os
import numpy as np
import matplotlib.pyplot as plt
from IO_Func import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=10,fontsize=16)

BaseDir = '../Numerical_Tests/1d_Standard_Soliton/En-Mass_Cp/'
f_list = Get_File_List(BaseDir)
StdGaussCollo = [fname for fname in f_list if fname.startswith("GL")]
PPClist = [fname for fname in f_list if fname.startswith("PPC")]
def GetEEvsDt(FoldNameList): 
    t_set = []
    ee_set = []
    for fname in FoldNameList:
        ffname = os.path.join(BaseDir,fname)
        err_file_name = Get_File_List(ffname)[0]
        if err_file_name.startswith('dt_'):
            data = np.load(os.path.join(ffname,err_file_name),allow_pickle=True).item()
            energyerr = np.array([singleengy-data['endt_energ'][0] for singleengy in data['endt_energ']])
            t_set.append(data['tau'])
            ee_set.append(max(abs(energyerr)))
    indice = np.argsort(np.array(t_set)).copy()
    t_set = np.sort(np.array(t_set))
    newarr = np.array(ee_set)
    newarr = newarr[indice]
    return t_set, newarr

t_set1, ee_set1 = GetEEvsDt(StdGaussCollo)
t_set2, ee_set2 = GetEEvsDt(PPClist)
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.loglog(t_set1,ee_set1,'o-',label='The proposed method')
ax.loglog(t_set2,ee_set2,'^--',label='The standard collcation method')
plt.grid(True,ls='--')
ax.set_xlabel('$\\tau$')
ax.set_ylabel('Energy error')
ax.legend(loc='upper left')
plt.show()

# SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
# fig.savefig(os.path.join(SavePath,'{}.pdf'.format('Example1EEvsdt')),dpi=600,format='pdf')
