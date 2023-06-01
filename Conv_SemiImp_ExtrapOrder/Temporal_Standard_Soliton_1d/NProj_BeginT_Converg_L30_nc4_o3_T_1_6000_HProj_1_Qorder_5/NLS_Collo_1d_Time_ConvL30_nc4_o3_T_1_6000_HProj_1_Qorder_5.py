from Collocation import CPNLS_SemiImpLin, Lap_Type
import os
import numpy as np
from Package_G1dMesh import Mesh1d
from Package_MyCode import FO
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
    'Temporal_Order2_Soliton': {
        'L': 30,
        '2': [60,70,80,90,100],
        '3': [20,25,30,35,40],
        '4': [8,12,14,16,20],
        'N': 1920,
        'T': 'np.pi/4'
    },
    'Temporal_Standard_Soliton': {
        'L': 30,
        '1': [60,70,80,90,100],
        '2': [60,70,80,90,100],
        '3': [20,25,30,35,40],
        '4': [16,24,48,64,80],
        'N': 6000, # spatial decomposition number
        'T': '1'
    },
    'Temporal_Bi_Soliton_T_Inv': {
        'L': 20,
        '2': [16,24,32,40,48,56],
        '3': [12,18,24,30,36,48],
        '4': [8,12,16,20,24,32],
        'N': 1024,
        'T': '2'
    },
    'Temporal_Bi_Soliton_T': {
        'L': 20,
        '1': [40,60],
        '2': [24,36,48,60,72],
        '3': [24,36,48,60,72],
        '4': [80,120,160,200],
        'N': 2048,
        'T': '2'
    },
    'ME_Bi_Soliton_T': {
        'L': 20,
        '3': [20],
        'N': 200,
        'T': '2'
    },
    'ME_Standard_Soliton': {
        'L': 40,
        '3': [20],
        'N': 400, # spatial decomposition number
        'T': '1'
    },
}
Lap_opt = Lap_Type.Diri_opt
# Mesh Param
Conv_type = 'Temporal_Standard_Soliton'
Md_type = 'BeginT'
dim = 1
L = N_sets_Dict[Conv_type]['L']
lp, rp = -L, L
order = 3
N_thres = 1e-12
ref_order = 1
N = N_sets_Dict[Conv_type]['N']
# Method Param
Tstr = N_sets_Dict[Conv_type]['T']
T = eval(Tstr)
quad_order = 5 # not used
x_save_ref = np.linspace(lp,rp,2**10+1)

for n_collo in [4]:
    save_lowest_index = 0
    print('Parameter of this example:')
    print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
    print('T = {}, N = {}'.format(Tstr, N))

    suffix = 'L{}_nc{}_o{}_T_{}_{}_HProj_{}_Qorder_{}'.format(L,n_collo,order,Tstr,N,ref_order,quad_order).replace('/',':')
    BaseDirPath = '/home/jiashhu/Collocation_NLS/Conv_SemiImp_ExtrapOrder/{}_{}d/NProj_{}_Converg_{}'.format(Conv_type,dim,Md_type,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/Semi_Implicit_Linear/',TargetDir=BaseDirPath,CodeName='NLS_Collo_1d_Time_Conv.py',suffix=suffix)

    N_Tsteps_sets = N_sets_Dict[Conv_type][str(n_collo)]

    for N_Tsteps in N_Tsteps_sets:
        if (save_lowest_index == 0) & (N_Tsteps == max(N_Tsteps_sets)):
            savefunc_obj = yxt_1d_out(x_save_ref,[T],[20],BaseDirPath)
            save_lowest_index += 1
        else:
            savefunc_obj = None
        dt = T/N_Tsteps
        # Newton iteration
        res_name = 'dt_{}_h_{}.npy'.format(N_Tsteps,N)

        # Generate 1d Periodic Mesh
        Mesh_Obj = Mesh1d(lp,rp,N,periodic=True)
        mesh     = Mesh(Mesh_Obj.ngmesh)
        myObj    = CPNLS_SemiImpLin(
            mesh, kappa = 2, p = 3,
            order = order, n_collocation = n_collo, 
            dt = dt, T = T, N_thres = N_thres, 
            ref_order = ref_order, quad_order = quad_order
        )
        myObj.Set_Lap_Type(Lap_opt)
        myObj.WeakForm4ImpLinExp()
        myObj.Sol_Setting(Exact_u_Dict['_'.join(Conv_type.split('_')[1:])])
        myObj.IniByProj()
        myObj.Solving(myObj.PP_Pic,so=savefunc_obj,Eig_opt=True)
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
        except:
            pass