import os
from Package_MyCode import FO
import numpy as np
import matplotlib.pyplot as plt

FoldName1 = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Temporal_Standard_Soliton_1d/Converg_L20_nc2_o3_T_1_1024_HProj_1_Qorder_5/dt_512_h_1024.npy'
FoldName2 = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Temporal_Standard_Soliton_1d/NoP_Converg_L20_nc2_o3_T_1_1024_HProj_1_Qorder_5/dt_512_h_1024.npy'
def GetErrvsDt(FoldName): 
    data = np.load(FoldName,allow_pickle=True).item()
    t_set = np.array(data['endt_T_set'])
    # intH1err = np.array(data['endt_H1err_ex'])
    intH1err = np.array(np.array(data['endt_energ'][1:])-data['endt_energ'][0])
    return t_set, np.abs(intH1err)

t_set1, intH1err1 = GetErrvsDt(FoldName1)
t_set2, intH1err2 = GetErrvsDt(FoldName2)
fig, ax = plt.subplots(1, 1, figsize=(7, 7))
ax.semilogy(t_set1,intH1err1,'ko',label = 'PP')
ax.semilogy(t_set2,intH1err2,'ro')
plt.legend()
plt.show()
