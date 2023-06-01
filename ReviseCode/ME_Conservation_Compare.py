import matplotlib.pyplot as plt
import numpy as np
from Package_MyCode import FO, PlotLineStyleObj
import matplotlib.gridspec as gridspec
import os

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=4)

BaseDir = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Spatial_Bi_Soliton_1d/'
CaseName = 'Converg_L20_nc2_o3_T_128_4096_HProj_3_Qorder_5_Lap_Dirichlet'
# CaseName = 'Converg_L20_nc2_o3_T_512_16384_HProj_3_Qorder_5_Lap_Dirichlet'
f_name = os.path.join(BaseDir, CaseName)
f_name1 = f_name
f_name2 = '/'.join(f_name.split('/')[:-1])+'/NoP_'+f_name.split('/')[-1]

# CaseName = 'Converg_L20_nc2_o3_T_128_8192_HProj_3_Qorder_5_Lap_Dirichlet'
# # CaseName = 'Converg_L20_nc2_o3_T_512_32768_HProj_3_Qorder_5_Lap_Dirichlet'
# f_name = os.path.join(BaseDir, CaseName)
# f_name3 = f_name
# f_name4 = '/'.join(f_name.split('/')[:-1])+'/NoP_'+f_name.split('/')[-1]

fig = plt.figure(figsize=(18,7))
gs = gridspec.GridSpec(19,1)

def GetMEDiff(f_name):
    dirs_list = [npyfile for npyfile in FO.Get_File_List(f_name) if npyfile.startswith('dt') and npyfile.endswith('_540.npy')]
    single_case = dirs_list[-1]
    # print('CaseName = {}'.format(single_case))
    # for single_case in dirs_list:
    #     if single_case.startswith('dt_'):
    npy_path = os.path.join(f_name,single_case)
    res = np.load(npy_path,allow_pickle=True).item()
    entH1err = np.append(res['endt_H1err_ex'][0],res['endt_H1err_ex'])
    entL2err = np.append(res['endt_L2err_ex'][0],res['endt_L2err_ex'])
    mass   = res['endt_massc']
    energy = res['endt_energ']
    tset   = np.append([0],res['endt_T_set'])
    mass_diff = [np.abs((mass[ii]-mass[0])) for ii in range(len(mass))]
    energy_diff = [np.abs((energy[ii]-energy[0])) for ii in range(len(energy))]

    return tset, mass_diff, energy_diff, entH1err, entL2err
        # print('energy is {}, mass is {}'.format(energy[0],mass[0]))


ax1 = fig.add_subplot(211)
ax1.set_position(gs[2:10].get_position(fig))

tset1, mass_diff1, energy_diff1, entH1err1, entL2err1 = GetMEDiff(f_name1)
tset2, mass_diff2, energy_diff2, entH1err2, entL2err2 = GetMEDiff(f_name2)
# tset3, mass_diff3, energy_diff3, entH1err3, entL2err3 = GetMEDiff(f_name3)
# tset4, mass_diff4, energy_diff4, entH1err4, entL2err4 = GetMEDiff(f_name4)
ax1.semilogy(tset1[1:],energy_diff1[1:],'-',label='Post-processing correction $\\tau = 2^{-5}$')
ax1.semilogy(tset1[1:],energy_diff2[1:],'--',label='Standard collocation $\\tau = 2^{-5}$')
# ax1.semilogy(tset3[1:],energy_diff3[1:],'-',label='Post-processing correction $\\tau = 2^{-6}$')
# ax1.semilogy(tset3[1:],energy_diff4[1:],'-.',label='Standard collocation $\\tau = 2^{-6}$')

ax2 = fig.add_subplot(212)
ax2.set_position(gs[11:].get_position(fig))

ax2.semilogy(tset1[1:],entH1err1[1:],'-',label='Energy Proposed')
ax2.semilogy(tset1[1:],entH1err2[1:],'--',label='Energy Standard')
# ax2.semilogy(tset3[1:],entH1err3[1:],'-',label='Energy Proposed2')
# ax2.semilogy(tset3[1:],entH1err4[1:],'-.',label='Energy Standard2')

ax1.set_ylim([1e-15,1e-4])
ax1.set_ylabel('Energy loss')
# ax1.set_xlabel('t')
ax2.set_xlabel('t')
ax2.set_ylabel('$H^1$ error')
ax1.legend(bbox_to_anchor=(0.5, 1.4), loc='upper center', ncol=2)
plt.show()

SavePath = '/Users/liubocheng/Documents/2022/Code-Collocation_NLS/ReviseCode/Pic'
fig.savefig(os.path.join(SavePath,'{}.pdf'.format('Example2EEvst2')),dpi=600,format='pdf')
