from Collocation import CPN_Converg_2d_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d, PeriodicSquare
from Package_MyCode import FO
from ngsolve import Mesh

# Method Param
kappa = -2
p = 3
n_collo = 3         # k=3 in Feng & Ma & Li
order = 2
N_thres = 1e-10
T = 0.1
N_Tsteps = 20
ref_order = 2
dim = 2
Conv_type = 'Spatial'

print('Parameter of this example:')
print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
print('T = {}, N_Tsteps = {}'.format(T, N_Tsteps))

suffix = 'nc{}_o{}_{}_{}_HProj_{}'.format(n_collo,order,N_Tsteps,N_thres,ref_order)
BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/{}/Converg{}d_{}'.format(Conv_type,dim,suffix)
if not os.path.exists(BaseDirPath):
    os.mkdir(BaseDirPath)

FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_2d_Spatial_Conv.py',suffix=suffix)

if order == 1:
    h_sets = [50,70,100]
elif order == 2:
    h_sets = [10,15,20,25,30]
elif order == 3:
    h_sets = [12,14,16,18,20]

for h_Spatial in h_sets:
    dt = T/N_Tsteps
    h = 1/h_Spatial
    # Newton iteration
    res_name = 'dt_{}_h_{}'.format(N_Tsteps,h_Spatial).replace('.','_')+'.npy'

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