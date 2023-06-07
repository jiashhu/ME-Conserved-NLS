import os
import sys
import json
import inspect
sys.path.append('/Users/liubocheng/Documents/2023/ME-Conserved-NLS/Packages')

import numpy as np
from Collocation import CPN_Converg_1d_Focus_Newt
from Package_G1dMesh import Mesh1d
from Package_MyCode import FO
from Package_ALE_Geometry import yxt_1d_out
from ngsolve import *
from Exact_Sol import *

def Main(dim,L,N_S,N_T,T,k_T,p_S,Output,Sol_Type,Ncore,Test_Name):
    '''
        dim: dimension of the problem
        L: half length of the truncated computational domain in 1d, [-L,L] 
        N_S: Number of intervals
        N_T: Number of time steps
        T: final evolution time
        k_T: Gauss collocation order
        p_S: spatial FEM order
        Output: decide whether print the numerical solution
            Optional: ["Err", "Sol", "All"]
        Sol_Type: decide which exact solution is used
    '''

    Lap_opt = 'Dirichlet'
    lp, rp = -L, L
    N_thres = 1e-10
    ref_order = 1
    dt = T/N_T

    print('Parameter of this example:')
    print('dim = {}, T_collo = {}, S_order = {}'.format(dim, k_T,p_S))
    print('L = {}, N_S = {}, T = {}, N_T = {}'.format(L, N_S, T, N_T))

    suffix = 'L{}_{}_nc{}_o{}_T_{}_{}'.format(L,N_S,k_T,p_S,T,N_T).replace('/',':')
    BaseDirPath = '../{}d_{}/{}/Conv_{}'.format(dim,Sol_Type,Test_Name,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    if Output in ["Err", "All"]:
        res_name = 'dt_{}_h_{}.npy'.format(N_T,N_S)

    if Output in ["Sol", "All"]:
        savefunc_obj = yxt_1d_out(np.linspace(lp,rp,2**10+1),[T],[20],BaseDirPath)
    else:
        savefunc_obj = None

    # Generate 1d Mesh
    Mesh_Obj = Mesh1d(lp,rp,N_S,periodic=False)
    mesh     = Mesh(Mesh_Obj.ngmesh)
    myObj    = CPN_Converg_1d_Focus_Newt(
        mesh, kappa = 2, p = 3,
        order = p_S, n_collocation = k_T, 
        dt = dt, T = T, N_thres = N_thres, 
        ref_order = ref_order
    )
    myObj.Set_Lap_Type(Lap_opt)
    myObj.WeakForm4Newton_Iter()
    myObj.Sol_Setting(Exact_u_Dict[Sol_Type])
    myObj.IniByProj()
    myObj.Solving(save_obj=savefunc_obj,ThreadNum=Ncore)
    Res_dict = {
        "endt_T_set":np.linspace(0,T,N_T+1)[1:],
        "endt_L2err_ex":myObj.endt_L2_ex_err_set,
        "endt_H1err_ex":myObj.endt_H1_ex_err_set,
        "int_L2err_ex":myObj.int_L2_ex_err_set,
        "int_H1err_ex":myObj.int_H1_ex_err_set,
        "endt_massc":myObj.endt_mas_c_set,
        "endt_energ":myObj.endt_eng_c_set,
        "int_massc":myObj.int_mas_c_set,
        "int_energ":myObj.int_eng_c_set,
        "int_T_set":myObj.int_tset,
        "endt_N_iter":myObj.endt_Newton_it,
        "endt_Param_iter":myObj.endt_alpbet_it,
        "endt_exactu_H1":myObj.endt_exactu_H1,
        "int_exactu_H1":myObj.int_exactu_H1,
    }
    try:
        print(max([max(myObj.endt_H1_ex_err_set),max(myObj.int_H1_ex_err_set)]))
    except:
        pass
    np.save(os.path.join(BaseDirPath,res_name),Res_dict,allow_pickle=True)

    try:
        xyt_Fig_Name = 'NumSol_Nspace_{}_NT_{}.npy'.format(N_S,N_T)
        np.save(os.path.join(BaseDirPath,xyt_Fig_Name),{
            'fval': savefunc_obj.Data_List,
            'tset': savefunc_obj.t_List,
            'xref': savefunc_obj.x_ref
        },allow_pickle=True)
    except:
        pass

if __name__ == "__main__":
    # Get the file path from the command-line argument
    file_path = sys.argv[1]

    # Read the JSON input from the file
    with open(file_path, "r") as json_file:
        data = json.load(json_file)

    # Get the parameter names from the Main function
    param_names = inspect.signature(Main).parameters.keys()
    matching_data = {k: v for k, v in data.items() if k in Main.__code__.co_varnames}
    data_dict = dict(data)
    for param_name in matching_data:
        # Get the values to test for the current parameter
        value = data[param_name]

        if not isinstance(value, (list, tuple)):
            data_dict[param_name] = value
        else:
            # Iterate over the values and call Main function
            for item in value:
                # Create a copy of the JSON data
                data_copy = dict(data_dict)
                
                # Update the current parameter value
                data_copy[param_name] = item
                
                # Call Main function with updated JSON data
                Main(**data_copy)

    Main(**matching_data)

