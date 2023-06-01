from Collocation import CPN_Converg_2d_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d, PeriodicSquare
from Package_MyCode import FO
from ngsolve import Mesh

# Method Param
dim = 2
kappa = -2
p = 3
n_collo = 4         # k=3 in Feng & Ma & Li
order = 3
N_thres = 1e-10
T = 0.1
ref_order = 1
N = 80

print('Parameter of this example:')
print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
print('T = {}, N = {}'.format(T, N))

suffix = 'nc{}_o{}_{}_{}_HProj_{}'.format(n_collo,order,N,N_thres,ref_order)
BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/Time/Converg{}d_{}'.format(dim,suffix)
if not os.path.exists(BaseDirPath):
    os.mkdir(BaseDirPath)

FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_2d_Time_Conv.py',suffix=suffix)

if n_collo == 2:
    N_Tsetps_sets = [46,48,50,52,54]
elif n_collo == 3:
    N_Tsetps_sets = [6,8,10,12,14]
elif n_collo == 4:
    N_Tsetps_sets = [3,4,5,6,7]

for N_Tsteps in N_Tsetps_sets:
    dt = T/N_Tsteps
    h = 1/N
    # Newton iteration
    res_name = 'dt_{}_h_{}'.format(N_Tsteps,N).replace('.','_')+'.npy'

    # Generate 1d Periodic Mesh
    Mesh_Obj = PeriodicSquare(h,periodic=True)
    mesh     = Mesh(Mesh_Obj.ngmesh)
    myObj    = CPN_Converg_2d_Newt(
        mesh, kappa = kappa, p = p, 
        order = order, n_collocation = n_collo, 
        dt = dt, T = T, N_thres = N_thres, 
        ref_order = ref_order, quad_order = 5
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