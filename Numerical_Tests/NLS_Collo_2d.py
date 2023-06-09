from Collocation import CPN_Converg_2d_Newt
import os
import sys
import json
import inspect
import numpy as np
from Package_G1dMesh import PeriodicSquareUniform
from IO_Func import Vtk_out
from ngsolve import *
from Exact_Sol import *


def Main(dim,L,N_S,N_T,T,k_T,p_S,Output,Sol_Type,Ncore,Test_Name,PPC=True):
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
        PPC: true for using post-processing correction
    '''

    Lap_opt = 'Dirichlet'
    lp, rp = -L, L
    N_thres = 1e-10
    ref_order = 2
    dt = T/N_T
    STD = not PPC

    print('Parameter of this example:')
    print('dim = {}, T_collo = {}, S_order = {}'.format(dim, k_T,p_S))
    print('L = {}, N_S = {}, T = {}, N_T = {}'.format(L, N_S, T, N_T))
    
    if PPC:
        prefix = "PPC"
    else:
        prefix = "GL"

    suffix = 'L{}_{}_nc{}_o{}_T_{}_{}'.format(L,N_S,k_T,p_S,T,N_T).replace('/',':')
    BaseDirPath = './{}d_{}/{}/{}_{}'.format(dim,Sol_Type,Test_Name,prefix,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    if Output in ["Err", "All"]:
        res_name = 'dt_{}_h_{}.npy'.format(N_T,N_S)

    if Output in ["Sol", "All"]:
        savefunc_obj = Vtk_out([T],[N_T],BaseDirPath)
    else:
        savefunc_obj = None
    
    Res_dict = {}
    # Generate 2d Mesh
    Mesh_Obj = PeriodicSquareUniform(N=N_S,L=L,periodic=False,Rot_opt=True)
    mesh     = Mesh(Mesh_Obj.ngmesh)
    # myObj    = CPN_Converg_2d_Newt(
    #     mesh, kappa = 2, p = 3, 
    #     order = p_S, n_collocation = k_T, 
    #     dt = dt, T = T, N_thres = N_thres, 
    #     ref_order = ref_order
    # )
    # myObj.Set_Lap_Type(Lap_opt)
    # myObj.WeakForm4Newton_Iter()
    # myObj.Sol_Setting(Exact_u_Dict[Sol_Type])
    # myObj.IniByProj()
    # myObj.Solving(save_obj=savefunc_obj,ThreadNum=Ncore,Eig_opt=STD)

    # Res_dict = {
    #     "endt_T_set":np.linspace(0,T,N_T+1)[1:],
    #     "endt_L2err_ex":myObj.endt_L2_ex_err_set,
    #     "endt_H1err_ex":myObj.endt_H1_ex_err_set,
    #     "int_L2err_ex":myObj.int_L2_ex_err_set,
    #     "int_H1err_ex":myObj.int_H1_ex_err_set,
    #     "endt_massc":myObj.endt_mas_c_set,
    #     "endt_energ":myObj.endt_eng_c_set,
    #     "int_massc":myObj.int_mas_c_set,
    #     "int_energ":myObj.int_eng_c_set,
    #     "int_T_set":myObj.int_tset,
    #     "endt_N_iter":myObj.endt_Newton_it,
    #     "endt_Param_iter":myObj.endt_alpbet_it,
    #     "endt_exactu_H1":myObj.endt_exactu_H1,
    #     "int_exactu_H1":myObj.int_exactu_H1,
    # }

    # if Output in ["Err", "All"]:
    #     np.save(os.path.join(BaseDirPath,res_name),Res_dict,allow_pickle=True)


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
    k = 0
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
                k = 1
    if k == 0:
        Main(**matching_data)