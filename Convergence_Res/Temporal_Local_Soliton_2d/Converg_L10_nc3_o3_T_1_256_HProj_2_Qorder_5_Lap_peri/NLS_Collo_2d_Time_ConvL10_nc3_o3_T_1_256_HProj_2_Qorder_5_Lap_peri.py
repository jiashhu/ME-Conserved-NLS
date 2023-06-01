from Collocation import CPN_Converg_2d_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d, PeriodicSquare, PeriodicSquareUniform
from Package_MyCode import FO
from ngsolve import *

Exact_u_Dict = {
    'Temporal_Local_Soliton': lambda t: 4/(2*exp(x+y-2*np.sqrt(2)*t)+exp(-(x+y-2*np.sqrt(2)*t)))*exp(1j*(t+(x+y)*np.sqrt(2)/2)),
    'Test': lambda t: 4/(2*exp(x+y-2*np.sqrt(2)*t)+exp(-(x+y-2*np.sqrt(2)*t)))*exp(1j*(t+(x+y)*np.sqrt(2)/2))
}
N_sets_Dict = {
    'Temporal_Local_Soliton': {
        'L': 10,
        '2': [9,15],
        '3': [2,3,4],
        '4': [2,3,4,5],
        'N': 256,
        'T': '1'
    },
}

# Method Param
dim = 2
Conv_type = 'Temporal_Local_Soliton'
kappa = 2
p = 3
n_collo = 2         # k=3 in Feng & Ma & Li
order = 3
N_thres = 1e-10
ref_order = 2
N_Spatial = N_sets_Dict[Conv_type]['N']
# Method Param
Tstr = N_sets_Dict[Conv_type]['T']
T = eval(Tstr)
L = N_sets_Dict[Conv_type]['L']
quad_order = 5 # not used
Lap_opt = 'peri'

for n_collo in [3]:
    print('Parameter of this example:')
    print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
    print('T = {}, N = {}'.format(Tstr, N_Spatial))

    suffix = 'L{}_nc{}_o{}_T_{}_{}_HProj_{}_Qorder_{}_Lap_{}'.format(L,n_collo,order,Tstr,N_Spatial,ref_order,quad_order,Lap_opt).replace('/',':')
    BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/{}_{}d/Converg_{}'.format(Conv_type,dim,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_2d_Time_Conv.py',suffix=suffix)

    N_Tsteps_sets = N_sets_Dict[Conv_type][str(n_collo)]

    for N_Tsteps in N_Tsteps_sets:
        dt = T/N_Tsteps
        h = 1/N_Spatial
        # Newton iteration
        res_name = 'dt_{}_h_{}'.format(N_Tsteps,N_Spatial).replace('.','_')+'.npy'

        # Generate 1d Periodic Mesh
        Mesh_Obj = PeriodicSquareUniform(N=N_Spatial,L=L,periodic=True,Rot_opt=True)
        mesh     = Mesh(Mesh_Obj.ngmesh)
        myObj    = CPN_Converg_2d_Newt(
            mesh, kappa = kappa, p = p, 
            order = order, n_collocation = n_collo, 
            dt = dt, T = T, N_thres = N_thres, 
            ref_order = ref_order, quad_order = 5
        )
        myObj.Set_Lap_Type(Lap_opt)
        myObj.WeakForm4Newton_Iter()
        myObj.Sol_Setting(Exact_u_Dict[Conv_type])
        myObj.IniByProj()
        myObj.Solving()
        Res_dict = {
            "endt_T_set":np.linspace(0,T,N_Tsteps+1)[1:],
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
            print(N_Spatial)
        except:
            pass
        np.save(os.path.join(BaseDirPath,res_name),Res_dict,allow_pickle=True)