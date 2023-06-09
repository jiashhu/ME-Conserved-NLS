import matplotlib.pyplot as plt
import numpy as np
from IO_Func import *
import matplotlib.gridspec as gridspec
import os

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle(markersize=4)

BaseDir = '../Numerical_Tests/1d_Bi_Soliton/LongTSol/'
CaseName = 'L20_640_nc2_o3_T_128_4096_0'
f_name1 = os.path.join(BaseDir,"GL_{}".format(CaseName))
f_name2 = os.path.join(BaseDir,"PPC_{}".format(CaseName))

fig = plt.figure(figsize=(18,7))
gs = gridspec.GridSpec(19,1)

def GetMEDiff(f_name):
    dirs_list = [npyfile for npyfile in Get_File_List(f_name) if npyfile.startswith('dt') and npyfile.endswith('_640.npy')]
    single_case = dirs_list[-1]
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


ax1 = fig.add_subplot(211)
ax1.set_position(gs[2:10].get_position(fig))

tset1, mass_diff1, energy_diff1, entH1err1, entL2err1 = GetMEDiff(f_name1)
tset2, mass_diff2, energy_diff2, entH1err2, entL2err2 = GetMEDiff(f_name2)
ax1.semilogy(tset1[1:],energy_diff1[1:],'-',label='Post-processing correction $\\tau = 2^{-5}$')
ax1.semilogy(tset1[1:],energy_diff2[1:],'--',label='Standard collocation $\\tau = 2^{-5}$')

ax2 = fig.add_subplot(212)
ax2.set_position(gs[11:].get_position(fig))

ax2.semilogy(tset1[1:],entH1err1[1:],'-',label='Energy Proposed')
ax2.semilogy(tset1[1:],entH1err2[1:],'--',label='Energy Standard')

ax1.set_ylim([1e-15,1e-4])
ax1.set_ylabel('Energy loss')
ax2.set_xlabel('t')
ax2.set_ylabel('$H^1$ error')
ax1.legend(bbox_to_anchor=(0.5, 1.4), loc='upper center', ncol=2)
plt.show()
