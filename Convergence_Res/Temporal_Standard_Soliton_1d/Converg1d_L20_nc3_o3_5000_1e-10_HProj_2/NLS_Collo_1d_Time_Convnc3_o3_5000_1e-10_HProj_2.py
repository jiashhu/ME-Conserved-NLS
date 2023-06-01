from Collocation import CPN_Converg_1d_Focus_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d
from Package_MyCode import FO
from ngsolve import Mesh

# Mesh Param
dim = 1
L = 20
lp, rp = -L, L
# Method Param
n_collo = 3         # k=3 in Feng & Ma & Li
order = 3
N_thres = 1e-10
ref_order = 2
T = 1
N = 5000

print('Parameter of this example:')
print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
print('T = {}, N = {}'.format(T, N))

suffix = 'nc{}_o{}_{}_{}_HProj_{}'.format(n_collo,order,N,N_thres,ref_order)
BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/Time/Converg{}d_{}'.format(dim,suffix)
if not os.path.exists(BaseDirPath):
    os.mkdir(BaseDirPath)

FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_1d_Time_Conv.py',suffix=suffix)

if n_collo == 2:
    N_Tsetps_sets = [60,70,80,90,100]
elif n_collo == 3:
    N_Tsetps_sets = [20,25,30,35,40]
elif n_collo == 4:
    N_Tsetps_sets = [8,12,14,16,20]

for N_Tsteps in N_Tsetps_sets:
    dt = T/N_Tsteps
    # Newton iteration
    res_name = 'dt_{}_h_{}.npy'.format(N_Tsteps,N)

    # Generate 1d Periodic Mesh
    Mesh_Obj = Mesh1d(lp,rp,N,periodic=True)
    mesh     = Mesh(Mesh_Obj.ngmesh)
    myObj    = CPN_Converg_1d_Focus_Newt(
        mesh, kappa = 2, p = 3,
        order = order, n_collocation = n_collo, 
        dt = dt, T = T, N_thres = N_thres, 
        ref_order = ref_order
    )
    # myObj.WeakForm4FixP_Iter()
    myObj.WeakForm4Newton_Iter()
    myObj.Sol_Setting()
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
        "int_T_set":myObj.int_tset
    }
    np.save(os.path.join(BaseDirPath,res_name),Res_dict,allow_pickle=True)