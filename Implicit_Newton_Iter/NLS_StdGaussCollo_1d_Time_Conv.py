'''
    The standard Gauss Collocation Method, Temporal convergence 
'''

from Collocation import CPN_Converg_1d_Focus_Newt, Lap_Type
import os
import numpy as np
from Package_G1dMesh import Mesh1d
from ngsolve import *
from Package_ALE_Geometry import yxt_1d_out

sech = lambda x_: 1/cosh(x_)
tanh = lambda x_: sinh(x_)/cosh(x_)
sol_M = 1.2
sol_N = 1
J = np.arctanh(2*sol_M*sol_N/(sol_M**2+sol_N**2))

chx = (exp(x)+exp(-x))/2
Exact_u_Dict = {
    'Order2_Soliton': lambda t: (4+3*(exp(8*1j*t)-1)/chx**2)/(
            2*chx+3/4*(cos(8*t)-1)/chx**3
        )*exp(1j*t),
    'Standard_Soliton': lambda t: 2/(exp(x+4*t)+exp(-x-4*t))*exp(-1j*(2*x+3*t)),
    'Bi_Soliton_T_Inv': lambda t: (
            exp(-1j*sol_M**2*(2-t))*sol_M*sech(sol_M*x) - exp(-1j*sol_N**2*(2-t))*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) + cos((sol_M**2-sol_N**2)*(2-t))*sech(sol_M*x)*sech(sol_N*x))
        ),
    'Bi_Soliton_T': lambda t: (
            exp(1j*sol_M**2*t)*sol_M*sech(sol_M*x) - exp(1j*sol_N**2*t)*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) + cos((sol_M**2-sol_N**2)*t)*sech(sol_M*x)*sech(sol_N*x))
        )
}
N_sets_Dict = {
    'Temporal_Standard_Soliton': {
        'L': 20,
        '2': np.array([8,32,128,512,2048]),
        '3': [20,25,30,35,40],
        '4': [8,12,14,16,20],
        'N': 1024, # spatial decomposition number
        'T': '1'
    }
}
Lap_opt = Lap_Type.Diri_opt
# Mesh Param
Conv_type = 'Temporal_Standard_Soliton'
dim = 1
L = N_sets_Dict[Conv_type]['L']
N = N_sets_Dict[Conv_type]['N']
Tstr = N_sets_Dict[Conv_type]['T']
lp, rp = -L, L
order = 3
N_thres = 1e-10
ref_order = 1
T = eval(Tstr)
quad_order = 5 # not used
x_save_ref = np.linspace(lp,rp,2**10+1)

for n_collo in [2]:
    save_lowest_index = 0
    print('Parameter of this example:')
    print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
    print('T = {}, N = {}'.format(Tstr, N))

    suffix = 'L{}_nc{}_o{}_T_{}_{}_HProj_{}_Qorder_{}'.format(L,n_collo,order,Tstr,N,ref_order,quad_order).replace('/',':')
    BaseDirPath = './Data_Convergence/{}_{}d/NoP_Converg_{}'.format(Conv_type,dim,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    N_Tsteps_sets = N_sets_Dict[Conv_type][str(n_collo)]

    for N_Tsteps in N_Tsteps_sets:
        # 如果从来没有存储过，并且时间步长是所有待选参数里面最小的那个
        if (save_lowest_index == 0) & (N_Tsteps == max(N_Tsteps_sets)):
            savefunc_obj = yxt_1d_out(x_save_ref,[T],[40],BaseDirPath)
            save_lowest_index += 1
        else:
            savefunc_obj = None
        dt = T/N_Tsteps
        # Newton iteration
        res_name = 'dt_{}_h_{}.npy'.format(N_Tsteps,N)

        # Generate 1d Periodic Mesh
        Mesh_Obj = Mesh1d(lp,rp,N,periodic=(Lap_opt=='peri'))
        mesh     = Mesh(Mesh_Obj.ngmesh)
        myObj    = CPN_Converg_1d_Focus_Newt(
            mesh, kappa = 2, p = 3,
            order = order, n_collocation = n_collo, 
            dt = dt, T = T, N_thres = N_thres, 
            ref_order = ref_order
        )
        myObj.Set_Lap_Type(Lap_opt)
        myObj.WeakForm4Newton_Iter()
        myObj.Sol_Setting(Exact_u_Dict['_'.join(Conv_type.split('_')[1:])])
        myObj.IniByProj()
        myObj.Solving(save_obj=savefunc_obj,Eig_opt=True)
        Res_dict = {
            "tau": dt,
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
            print(N_Tsteps)
        except:
            pass
        np.save(os.path.join(BaseDirPath,res_name),Res_dict,allow_pickle=True)

        try:
            xyt_Fig_Name = 'NumSol_Nspace_{}_NT_{}.npy'.format(N,N_Tsteps)
            np.save(os.path.join(BaseDirPath,xyt_Fig_Name),{
                'fval': savefunc_obj.Data_List,
                'tset': savefunc_obj.t_List,
                'xref': savefunc_obj.x_ref
            },allow_pickle=True)
            print('NumSol saved Successful!')
        except:
            pass