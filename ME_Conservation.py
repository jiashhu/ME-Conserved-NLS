import matplotlib.pyplot as plt
import numpy as np
from Package_MyCode import FO
import os

f_name = 'Convergence_Res/ME_Conserv'
dirs_list = FO.Get_File_List(f_name)
for single_case in dirs_list:
    data_fname = FO.Get_File_List(os.path.join(f_name,single_case))
    assert len(data_fname)==1
    npy_path = os.path.join(f_name,single_case,data_fname[0])
    res = np.load(npy_path,allow_pickle=True).item()
    mass = res['endt_massc']
    energy = res['endt_energ']
    tset = res['endt_T_set']

    mass_diff = [np.real(mass[ii]-mass[0]) for ii in range(len(mass))]
    energy_diff = [np.real(energy[ii]-energy[0]) for ii in range(len(energy))]

    fig, ax = plt.subplots(1, 2, figsize=(12, 9))
    ax[0].plot(tset,mass_diff,'-')
    ax[1].plot(tset,energy_diff,'-')
    plt.show()

