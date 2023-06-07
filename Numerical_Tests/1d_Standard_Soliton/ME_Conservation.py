import matplotlib.pyplot as plt
import numpy as np
from Package_MyCode import FO, PlotLineStyleObj
import os

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=4)

f_name = '/Users/hjs/Documents/2022/Collocation_NLS/Convergence_Res/ME_Conserv/ME_Standard_Soliton_1d/Converg_L40_nc3_o3_T_1_400_HProj_1_Qorder_5'
fig, ax = plt.subplots(2, 1, figsize=(7, 7))
dirs_list = FO.Get_File_List(f_name)
for single_case in dirs_list:
    if single_case.startswith('dt_20'):
        npy_path = os.path.join(f_name,single_case)
        res = np.load(npy_path,allow_pickle=True).item()
        mass   = res['endt_massc']
        energy = res['endt_energ']
        tset   = np.append([0],res['endt_T_set'])

        mass_diff = [np.real(mass[ii]-mass[0]) for ii in range(len(mass))]
        energy_diff = [np.real(energy[ii]-energy[0]) for ii in range(len(energy))]
        ax[0].plot(tset,mass_diff,'o-',label='$M_h(t)-M_h(0)$')
        ax[0].plot(tset,energy_diff,'>--',label='$E_h(t)-E_h(0)$')

        N_iter = np.append([np.nan],res['endt_N_iter'])
        P_iter = np.append([np.nan],res['endt_Param_iter'])
        ax[1].plot(tset,N_iter,'o-',label='Collocation FEM Iteration')
        ax[1].plot(tset,P_iter,'>--',label='Parameter Iteration')

ax[0].legend(loc='upper right')
ax[0].set_ylim([-1e-14,1e-14])
ax[1].legend()
ax[1].set_xlabel('t')
ax[1].set_xlim([0,1])
ax[0].set_xlim([0,1])
# plt.show()

SavePath = '/Users/hjs/Documents/2022/Collocation_NLS/Draw_Pic/Pic'
fig.savefig(os.path.join(SavePath,'{}.eps'.format('ME_ITER_Std_Soliton_T_1d')),dpi=600,format='eps')
