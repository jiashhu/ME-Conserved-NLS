import matplotlib.pyplot as plt
import numpy as np
from Package_MyCode import FO, PlotLineStyleObj
import os

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=4)

# Proj = 'NoP'
Proj = ''
f_name = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_T_Inv_1d/NoP_Converg_L20_nc3_o2_T_256_8192_HProj_3_Qorder_5_Lap_Dirichlet'
f_name2 = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_T_Inv_1d/Converg_L20_nc3_o2_T_256_8192_HProj_3_Qorder_5_Lap_Dirichlet'
fig, ax = plt.subplots(2, 1, figsize=(7, 7))
dirs_list = FO.Get_File_List(f_name)
dirs_list2 = FO.Get_File_List(f_name2)
single_case = [vitem for vitem in dirs_list if vitem.startswith('dt_')]
single_case2 = [vitem for vitem in dirs_list2 if vitem.startswith('dt_')]

npy_path1 = os.path.join(f_name,single_case[0])
npy_path2 = os.path.join(f_name2,single_case2[0])
index = 0

for npy_path in [npy_path1, npy_path2]:
    res = np.load(npy_path,allow_pickle=True).item()
    entH1err = np.append(res['endt_H1err_ex'][0],res['endt_H1err_ex'])
    entL2err = np.append(res['endt_L2err_ex'][0],res['endt_L2err_ex'])
    mass   = res['endt_massc']
    energy = res['endt_energ']
    tset   = np.append([0],res['endt_T_set'])

    mass_diff = [np.real((mass[ii]-mass[0])) for ii in range(len(mass))]
    energy_diff = [np.real((energy[ii]-energy[0])) for ii in range(len(energy))]
    print('energy is {}, mass is {}'.format(energy[0],mass[0]))
    ax[index].plot(tset,mass_diff,'-',label='$M_h(t)-M_h(0)$')
    ax[index].plot(tset,energy_diff,'--',label='$E_h(t)-E_h(0)$')
    ax[index].legend(loc='upper right')
    index = index + 1
    
ax[1].set_xlabel('t')
plt.show()

    # SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
    # fig.savefig(os.path.join(SavePath,'{}.pdf'.format('{}_Bi_Soliton_T_1d_Inv'.format(Proj))),dpi=600,format='pdf')
