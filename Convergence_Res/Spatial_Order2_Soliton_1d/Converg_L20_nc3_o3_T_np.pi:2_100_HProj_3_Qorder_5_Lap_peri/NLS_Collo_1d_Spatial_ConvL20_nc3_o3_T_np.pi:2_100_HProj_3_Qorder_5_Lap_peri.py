from Collocation import CPN_Converg_1d_Focus_Newt, Lap_Type
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
Lap_opt = Lap_Type.peri_opt

chx = (exp(x)+exp(-x))/2
Exact_u_Dict = {
    'Spatial_Order2_Soliton': lambda t: (4+3*(exp(8*1j*t)-1)/chx**2)/(
            2*chx+3/4*(cos(8*t)-1)/chx**3
        )*exp(1j*t),
    'Spatial_Standard_Soliton': lambda t: 2/(exp(x+4*t)+exp(-x-4*t))*exp(-1j*(2*x+3*t)),
    'Spatial_Bi_Soliton': lambda t: (
            exp(1j*sol_M**2*t)*sol_M*sech(sol_M*x) - exp(1j*sol_N**2*t)*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) + cos((sol_M**2-sol_N**2)*t)*sech(sol_M*x)*sech(sol_N*x))
        ),
    'Spatial_Bi_Soliton_T_Inv': lambda t: (
            exp(-1j*sol_M**2*(2-t))*sol_M*sech(sol_M*x) - exp(-1j*sol_N**2*(2-t))*sol_N*sech(sol_N*x))/(
            np.cosh(J) - np.sinh(J)*(tanh(sol_M*x)*tanh(sol_N*x) + cos((sol_M**2-sol_N**2)*(2-t))*sech(sol_M*x)*sech(sol_N*x))
        )
}
N_sets_Dict = {
    'Spatial_Order2_Soliton': {
        'L': 20,
        '1': [480, 720, 960, 1200, 1440],
        '2': [240, 360, 480, 600, 720],
        '3': [120, 180, 240, 300, 360],
        '4': [1920],
        'Tstr': 'np.pi/2',
        'N_T': 100
    },
    'Spatial_Standard_Soliton': {
        'L': 20,
        '1': [480, 720, 960, 1200, 1440],
        '2': [240, 360, 480, 600, 720],
        '3': [120, 180, 240, 300, 360],
        'Tstr': '1',
        'N_T': 1000
    },
    'Spatial_Bi_Soliton': {
        'L': 20,
        '1': [480, 720, 960, 1200, 1440],
        '2': [240, 360, 480, 600, 720],
        '3': [120, 180, 240, 300, 360],
        'Tstr': '2',
        'N_T': 100
    },
    'Spatial_Bi_Soliton_T_Inv': {
        'L': 20,
        '1': [480, 720, 960, 1200, 1440],
        '2': [240, 360, 480, 600, 720],
        '3': [120, 180, 240, 300, 360],
        'Tstr': '2',
        'N_T': 100
    }
}
# unchanged parameters
dim = 1
Conv_type = 'Spatial_Order2_Soliton'
n_collo = 3         # k=3 in Feng & Ma & Li
N_thres = 1e-10
kappa = 2
p = 3               # nonlinear: kappa*|u|^(p-1)
ref_order = 3
quad_order = 5
L = N_sets_Dict[Conv_type]['L']
Tstr =  N_sets_Dict[Conv_type]['Tstr']
N_Tsteps = N_sets_Dict[Conv_type]['N_T']
# unchanged parameters
lp, rp = -L, L
T = eval(Tstr)
x_save_ref = np.linspace(lp,rp,2**10+1)

for order in [3]:
    save_lowest_index = 0
    # Mesh Param
    print('Parameter of this example:')
    print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}, Laplace option is {}'.format(dim, n_collo,order, N_thres, ref_order, Lap_opt))
    print('T = {}, N_Tsteps = {}'.format(Tstr, N_Tsteps))

    suffix = 'L{}_nc{}_o{}_T_{}_{}_HProj_{}_Qorder_{}_Lap_{}'.format(L,n_collo,order,Tstr,N_Tsteps,ref_order,quad_order,Lap_opt).replace('/',':')
    BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/{}_{}d/Converg_{}'.format(Conv_type,dim,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_1d_Spatial_Conv.py',suffix=suffix)

    N_sets = N_sets_Dict[Conv_type][str(order)]

    for N_Spatial in N_sets:
        if (save_lowest_index == 0) & (N_Spatial == max(N_sets)):
            savefunc_obj = yxt_1d_out(x_save_ref,[T],[20],BaseDirPath)
            save_lowest_index += 1
        else:
            savefunc_obj = None
        dt = T/N_Tsteps
        N = N_Spatial
        # Newton iteration
        res_name = 'dt_{}_h_{}.npy'.format(N_Tsteps,N)

        # Generate 1d Periodic Mesh
        Mesh_Obj = Mesh1d(lp,rp,N,periodic=True)
        mesh     = Mesh(Mesh_Obj.ngmesh)
        myObj    = CPN_Converg_1d_Focus_Newt(
            mesh, kappa = kappa, p = p, order = order, 
            n_collocation = n_collo, 
            dt = dt, T = T, N_thres = N_thres, quad_order = quad_order
        )
        myObj.Set_Lap_Type(Lap_opt)
        myObj.WeakForm4Newton_Iter()
        myObj.Sol_Setting(Exact_u_Dict[Conv_type])
        myObj.IniByProj()
        myObj.Solving(save_obj=savefunc_obj)

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
        try:
            xyt_Fig_Name = 'NumSol_Nspace_{}_NT_{}.npy'.format(N_Spatial,N_Tsteps)
            np.save(os.path.join(BaseDirPath,xyt_Fig_Name),{
                'fval': savefunc_obj.Data_List,
                'tset': savefunc_obj.t_List,
                'xref': savefunc_obj.x_ref
            },allow_pickle=True)
        except:
            pass