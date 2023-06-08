import numpy as np
import os
from IO_Func import *
import json
import sys

def GetData(npyfile):
    InfH1_err = 0
    res = np.load(os.path.join(npyfile),allow_pickle=True).item()
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
            InfH1_err = endt_H1err
            max_t_index = endt_max_index
            max_t = res['endt_T_set'][max_t_index]
        else:
            InfH1_err = int_H1err
            max_t_index = int_max_index
            max_t = res['int_T_set'][max_t_index]
    except:
        pass

    Newton_its = None
    Params_its = None
    try:
        Newton_its = res['endt_N_iter']
        Params_its = res['endt_Param_iter']
        print('max Newton iters is {}, max Param iters is {}'.format(
            max(Newton_its), max(Params_its)
        ))  
    except:
        pass
    
    res_back = {
        "InfH1_err": np.real(InfH1_err),
        "mass": res['endt_massc'],
        "energy": res['endt_energ'], 
        "N_its": Newton_its, 
        "P_its": Params_its
    }
    return res_back

def Gen_Name(L,N_S,k_T,p_S,T,N_T,PPC):
    if PPC:
        prefix = "PPC"
    else:
        prefix = "GL"

    suffix = 'L{}_{}_nc{}_o{}_T_{}_{}'.format(L,N_S,k_T,p_S,T,N_T).replace('/',':')
    Fname = '{}_{}'.format(prefix,suffix)
    return Fname

if __name__ == "__main__":
    # Get the file path from the command-line argument
    file_path = sys.argv[1]

    Param_Names = [fname for fname in Get_File_List(file_path) 
                   if fname.endswith('.json') and not fname.startswith('Sum')]
    
    Sum_Dict = {
        "Converg": {}
    }

    for Param in sorted(Param_Names):
        # Read the JSON input from the file
        with open(os.path.join(file_path,Param), "r") as json_file:
            data = json.load(json_file)
        ErrDict = {
            "param": [],
            "err": []
        }
        varkey = [vkey for vkey, value in data.items() if isinstance(value,(list,tuple))][0]
        for value in data[varkey]:
            ErrDict["param"].append(value)
            dict_copy = dict(data)
            dict_copy[varkey] = value
            matching_data = {k: v for k, v in dict_copy.items() if k in Gen_Name.__code__.co_varnames}
            folder = os.path.join(file_path,Gen_Name(**matching_data))
            errs_name = [fname for fname in Get_File_List(folder) 
                         if fname.startswith('dt')]
            res = GetData(os.path.join(folder,errs_name[0]))
            ErrDict["err"].append(res["InfH1_err"])

        Sum_Dict["Converg"][Param.split('.')[0]] = ErrDict

    # Specify the file path
    f_path = "Sum_{}.json".format(file_path.split('/')[-1])
    save_path = os.path.join(file_path,f_path)

    # Write the dictionary to the JSON file
    with open(save_path, "w") as json_file:
        json.dump(Sum_Dict, json_file, indent=4)