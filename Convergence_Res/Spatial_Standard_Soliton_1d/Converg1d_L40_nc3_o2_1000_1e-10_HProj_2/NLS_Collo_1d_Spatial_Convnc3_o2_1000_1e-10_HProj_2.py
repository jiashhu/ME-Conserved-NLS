from Collocation import CPN_Converg_1d_Focus_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d
from Package_MyCode import FO
from ngsolve import Mesh

# Mesh Param
kappa = 2
p = 3
L = 40
lp, rp = -L, L
# Method Param
n_collo = 3         # k=3 in Feng & Ma & Li
order = 2
N_thres = 1e-10
T = 1
N_Tsteps = 1000
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

FO.Copy_Code(OriginDir='/home/jiashhu/Collocation_NLS/',TargetDir=BaseDirPath,CodeName='NLS_Collo_1d_Spatial_Conv.py',suffix=suffix)

if order == 1:
    N_sets = [800,1200,1600,2000,2400]
elif order == 2:
    N_sets = [120,180,240,300,360]
elif order == 3:
    N_sets = [60,90,120,150,180]

for N_Spatial in N_sets:
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
        dt = dt, T = T, N_thres = N_thres
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