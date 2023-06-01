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
# f_name = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/NoP_Converg_L30_nc3_o3_T_512_4096_HProj_3_Qorder_5_Lap_Dirichlet'
# f_name = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_T_Inv_1d/Converg_L20_nc3_o2_T_32_1024_HProj_3_Qorder_5_Lap_Dirichlet'
# f_name = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/Converg_L30_nc3_o3_T_512_4096_HProj_3_Qorder_5_Lap_Dirichlet'
# f_name = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/Converg_L20_nc3_o3_T_512_4096_HProj_3_Qorder_5_Lap_Dirichlet'
BaseDir = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/'
CaseName = 'NoP_Converg_L20_nc2_o3_T_512_16384_HProj_3_Qorder_5_Lap_Dirichlet'
f_name = os.path.join(BaseDir,CaseName)
fig, ax = plt.subplots(2, 1, figsize=(7, 7))
dirs_list = FO.Get_File_List(f_name)
for single_case in dirs_list:
    if single_case.startswith('dt_'):
        npy_path = os.path.join(f_name,single_case)
        res = np.load(npy_path,allow_pickle=True).item()
        entH1err = np.append(res['endt_H1err_ex'][0],res['endt_H1err_ex'])
        entL2err = np.append(res['endt_L2err_ex'][0],res['endt_L2err_ex'])
        mass   = res['endt_massc']
        energy = res['endt_energ']
        tset   = np.append([0],res['endt_T_set'])

        mass_diff = [np.real((mass[ii]-mass[0])) for ii in range(len(mass))]
        energy_diff = [np.real((energy[ii]-energy[0])) for ii in range(len(energy))]
        print('energy is {}, mass is {}'.format(energy[0],mass[0]))
        ax[0].plot(tset,mass_diff,'-',label='$M_h(t)-M_h(0)$')
        ax[0].plot(tset,energy_diff,'--',label='$E_h(t)-E_h(0)$')

        ax[1].plot(tset,entH1err,'-',label='H1err')
        ax[1].plot(tset,entL2err,'-',label='L2err')

ax[0].legend(loc='upper right')
# ax[0].set_ylim([-1e-13,1e-13])
ax[1].legend()
ax[1].set_xlabel('t')
ax[0].set_xlim([0,128])
# ax[1].set_xlim([0,128])
# ax[1].set_ylim([0,1])
plt.show()

# SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
# fig.savefig(os.path.join(SavePath,'{}.pdf'.format('{}_Bi_Soliton_T_1d_Inv'.format(Proj))),dpi=600,format='pdf')
