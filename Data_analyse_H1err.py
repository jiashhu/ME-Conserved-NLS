import numpy as np
import os
from Package_MyCode import FO

BaseDir_path = '/home/jiashhu/Collocation_NLS/Convergence_Res/Spatial_Local_Soliton_2d/Converg_L10_nc3_o3_T_1_256_HProj_2_Qorder_5'
npy_list = [fname for fname in FO.Get_File_List(BaseDir_path) if fname.endswith('npy') & fname.startswith('dt')]

EndtErr = []
for npyfile in npy_list:
    print('filename is {}'.format(npyfile))
    res = np.load(os.path.join(BaseDir_path,npyfile),allow_pickle=True).item()
    NT_endt_steps = len(res['endt_T_set'])
    NT_int_steps = len(res['int_T_set'])
    endt_H1err = max((res['endt_H1err_ex'][:NT_endt_steps]))
    endt_max_index = np.argmax((res['endt_H1err_ex'][:NT_endt_steps]))
    int_H1err = max((res["int_H1err_ex"][:NT_int_steps]))
    int_max_index = np.argmax((res["int_H1err_ex"][:NT_int_steps]))
    
    endt_L2err = max((res['endt_L2err_ex'][:NT_endt_steps]))
    int_L2err = max((res["int_L2err_ex"][:NT_int_steps]))

    try:
        if endt_H1err>int_H1err:
            H1_err = endt_H1err
            max_t_index = endt_max_index
            max_t = res['endt_T_set'][max_t_index]
        else:
            H1_err = int_H1err
            max_t_index = int_max_index
            max_t = res['int_T_set'][max_t_index]
    except:
        pass

    print("H1err is {} at t = {}, L2err is {}".format(H1_err,max_t,max(endt_L2err,int_L2err)))

    EndtErr.append(H1_err)
    mass_np = res['endt_massc']
    mass_diff = np.std(mass_np)
    mass_max_diff = np.max(abs(mass_np-np.mean(mass_np)))
    energ_np = res['endt_energ']
    energ_diff = np.std(energ_np)
    energ_max_diff = np.max(abs(energ_np-np.mean(energ_np)))
    print('endtH1err is {}, H1err is {}, mass std is {}, energy std is {}, max_diffs are {} and {}'.format(
        endt_H1err, H1_err, mass_diff, energ_diff, mass_max_diff, energ_max_diff
    ))
    try:
        endt_exactu_H1 = np.array(res['endt_exactu_H1'][:NT_endt_steps])
        int_exactu_H1 = np.array(res['int_exactu_H1'][:NT_int_steps])
        endt_relative_H1 = max(np.array(res['endt_H1err_ex'][:NT_endt_steps])/endt_exactu_H1)
        int_relative_H1 = max(np.array(res['int_H1err_ex'][:NT_int_steps])/int_exactu_H1)
        H1_rel_err = max(endt_relative_H1,int_relative_H1)
        tmp = np.argmax(endt_exactu_H1)
        print('Endt max exactu H1 {} at index {} with t = {}'.format(max(endt_exactu_H1),tmp,res['endt_T_set'][tmp]))
        print('H1err relative error is {}'.format(H1_rel_err))
    except:
        pass

    try:
        Newton_its = res['endt_N_iter']
        Params_its = res['endt_Param_iter']
        print('max Newton iters is {}, max Param iters is {}'.format(
            max(Newton_its), max(Params_its)
        ))  
    except:
        pass

[print(err.real) for err in EndtErr]