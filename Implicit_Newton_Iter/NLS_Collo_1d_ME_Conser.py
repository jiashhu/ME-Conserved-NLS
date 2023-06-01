from Collocation import CPN_Converg_1d_Focus_Newt
import os
import numpy as np
from Package_G1dMesh import Mesh1d
from Package_MyCode import FO
from ngsolve import Mesh

# for all time collocation degree, using fem p=3 for tau=h=0.2 to T = 2
for index_012 in range(3):
    # Mesh Param
    kappa = 2
    p = 3
    L = 40
    lp, rp = -L, L
    # Method Param
    n_collo = index_012+2       # k = 2,3,4
    order = 3
    N_thres = 1e-10
    T = 2
    N_Tsteps = 10
    ref_order = 1
    dim = 1
    Expr_type = 'ME_Conserv'

    print('Parameter of this example:')
    print('dim = {}, n_collo = {}, order = {}, N_thres = {}, ref_order = {}'.format(dim, n_collo,order, N_thres, ref_order))
    print('T = {}, N_Tsteps = {}'.format(T, N_Tsteps))

    suffix = 'L{}_nc{}_o{}_T_{}_{}_{}_HProj_{}'.format(L,n_collo,order,T,N_Tsteps,N_thres,ref_order)
    BaseDirPath = '/home/jiashhu/Collocation_NLS/Convergence_Res/{}/Converg{}d_{}'.format(Expr_type,dim,suffix)
    if not os.path.exists(BaseDirPath):
        os.makedirs(BaseDirPath)

    N_sets = [400]    # 2*L/0.2

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
            "int_T_set":myObj.int_tset,
            "endt_N_iter":myObj.endt_Newton_it,
            "endt_Param_iter":myObj.endt_alpbet_it
        }
        np.save(os.path.join(BaseDirPath,res_name),Res_dict,allow_pickle=True)