from Collocation import CPN_Converg_2d_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d, PeriodicSquare, PeriodicSquareUniform
from Package_MyCode import FO
from ngsolve import *

Exact_u_Dict = {
    'Spatial_Local_Soliton': lambda t: 4/(2*exp(x+y-2*np.sqrt(2)*t)+exp(-(x+y-2*np.sqrt(2)*t)))*exp(1j*(t+(x+y)*np.sqrt(2)/2)),
}
N_sets_Dict = {
    'Spatial_Local_Soliton': {
        'L': 10,
        '1': [64, 96, 128,160,192],
        '2': [32, 48, 64, 80, 96],
        '3': [16, 24, 32, 40, 48]
    },
}
# unchanged parameters
dim = 2
Conv_type = 'Spatial_Local_Soliton'
n_collo = 3         # k=3 in Feng & Ma & Li
N_thres = 1e-10
kappa = 2
p = 3               # nonlinear: kappa*|u|^(p-1)
ref_order = 2
quad_order = 5
L = N_sets_Dict[Conv_type]['L']
Tstr = '0.5'
N_Tsteps = 40
# unchanged parameters
T = eval(Tstr)
Lap_opt = 'peri'

for order in [1,2]:
    # Mesh Param
    print('Parameter of this example:')
    print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
    print('T = {}, N_Tsteps = {}'.format(Tstr, N_Tsteps))

    suffix = 'L{}_nc{}_o{}_T_{}_{}_HProj_{}_Qorder_{}'.format(L,n_collo,order,Tstr,N_Tsteps,ref_order,quad_order).replace('/',':')
    BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/{}_{}d/Converg_{}'.format(Conv_type,dim,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_2d_Spatial_Conv.py',suffix=suffix)

    N_sets = N_sets_Dict[Conv_type][str(order)]

    for N_Spatial in N_sets:
        dt = T/N_Tsteps
        # Newton iteration
        res_name = 'dt_{}_NSpatial_{}'.format(N_Tsteps,N_Spatial).replace('.','_')+'.npy'

        # Generate 1d Periodic Mesh
        Mesh_Obj = PeriodicSquareUniform(N=N_Spatial,L=L,periodic=True,Rot_opt=True)
        mesh     = Mesh(Mesh_Obj.ngmesh)
        myObj    = CPN_Converg_2d_Newt(
            mesh, kappa = kappa, p = p, 
            order = order, n_collocation = n_collo, 
            dt = dt, T = T, N_thres = N_thres, 
            ref_order = ref_order, quad_order = quad_order
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