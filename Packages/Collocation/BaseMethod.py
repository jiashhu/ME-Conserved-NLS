from ngsolve import *
from Collocation import GL
import numpy as np
from typing import Callable
from IO_Func import Vtk_out, LogTime, Printer
from ngsolve.krylovspace import CGSolver
import time

class GatherData():
    def __init__(self,IO) -> None:
        self.iter = 0
        self.residue = 0
        self.IO = IO

    def collect(self,iters,residue):
        self.iter = max(self.iter,iters)
        self.residue = max(self.residue,residue)
        if self.IO:
            print("{} CG ends with iters {} and residue {}".format(LogTime(),iters,residue))
    
    def print(self):
        print('max CG iters = {} and max residue = {}'.format(self.iter,self.residue))

MatInvOpt_List = ['Iter','Direct']
MatInvOpt = 'Iter'

class Coll_Proj_NLS_FixP():
    def __init__(self,mesh,kappa,p,order,n_collocation,dt,T,thres,quad_order,ref_order) -> None:
        self.mesh = mesh
        self.order = order
        self.dt = dt       
        self.n_collo = n_collocation
        self.Quad_order = quad_order
        self.NLS_Threshold = thres
        self.Param_Tolerance = thres
        self.Param_iters = 0
        self.kappa = kappa
        self.p = p
        self.T = T
        self.t = 0

        self.fesc_peri = Periodic(H1(self.mesh, order=self.order, complex=True))
        self.fesc_peri_ref  = Periodic(H1(self.mesh, order=self.order+ref_order, complex=True))
        self.fesc_peri_Lag = self.fesc_peri*NumberSpace(self.mesh, complex=True)   # periodic cannot avoid constant
        self.fesc_Dirichlet = H1(self.mesh, order=self.order, complex=True, dirichlet=self.mesh.Boundaries('.*'))
        self.fesc_Diri_ref = H1(self.mesh, order=self.order+ref_order, complex=True, dirichlet=self.mesh.Boundaries('.*'))
        
    def Set_Lap_Type(self,Lap_opt):
        self.Lap_opt = Lap_opt
        if self.Lap_opt == 'peri':
            self.fesc_working = self.fesc_peri
            self.fesc_working_ref = self.fesc_peri_ref
            # for solving Laplace equation
            self.GeneInvM_Peri_Lap_LagMul()
        elif self.Lap_opt == 'Dirichlet':
            self.fesc_working = self.fesc_Dirichlet
            self.fesc_working_ref = self.fesc_Diri_ref
            self.GeneInv_Lap()

        # every collocation point is assigned with a spatial fem func
        self.fes_ALL = self.fesc_working**self.n_collo
        self.U_N = self.fes_ALL.TrialFunction()
        self.V_N = self.fes_ALL.TestFunction()
        # For Newton iterations
        self.U_N_old = GridFunction(self.fes_ALL)
        self.U_N_update = GridFunction(self.fes_ALL)
        self.lhs_Niter = BilinearForm(self.fes_ALL)
        self.rhs_Niter = LinearForm(self.fes_ALL)
        # tn or updated fem sol at tn+1
        self.un = GridFunction(self.fesc_working)
        self.un_aux = GridFunction(self.fesc_working)
        self.eu_ref = GridFunction(self.fesc_working_ref)

        # GL_differential_matrix 
        # 初值对应  与  配点对应  与  多项式插值系数
        self.D_mat, self.D_row, self.extrap_coeff, self.collo_inner_nodes = self.GenerateD_mat()
        # FixPoint iterations
        self.NLS_iters = 0
        # Nonlinear terms
        self.F_norm2 = lambda norm: 2*self.kappa/(self.p+1)*norm**(self.p+1)
        self.f_norm2 = lambda norm: self.kappa*norm**(self.p-1)
        self.df_norm2 = lambda norm: (self.p-1)/2*self.kappa*norm**(self.p-3)
        self.alpha, self.beta = 1,0
        self.IntOmeg = lambda expr: Integrate(expr,self.mesh,element_wise=False,order=self.Quad_order)
        self.L2Omeg = lambda expr: self.IntOmeg(InnerProduct(expr,expr))
        self.H1semi = lambda expr: self.IntOmeg(InnerProduct(grad(expr),grad(expr)))
        self.H1Omeg = lambda expr: self.IntOmeg(InnerProduct(expr,expr) + InnerProduct(grad(expr),grad(expr)))

    def GeneInv_Lap(self):
        self.Lap_Diri_lhs = BilinearForm(self.fesc_Dirichlet)
        self.Lapu, self.Lapv = self.fesc_Dirichlet.TnT()
        self.Lap_Diri_lhs += -InnerProduct(grad(self.Lapu),grad(self.Lapv))*dx
        self.Lap_Diri_inv = self.Lap_Diri_lhs.Assemble().mat.Inverse(self.fesc_Dirichlet.FreeDofs(),inverse='pardiso')

    def GeneInvM_Peri_Lap_LagMul(self):
        self.Lap_Lag_lhs = BilinearForm(self.fesc_peri_Lag)
        self.Lapu, self.LagMul = self.fesc_peri_Lag.TrialFunction()
        self.Lapv, self.LagTest = self.fesc_peri_Lag.TestFunction()
        self.Lap_Lag_lhs += -InnerProduct(grad(self.Lapu),grad(self.Lapv))*dx + self.LagMul*self.Lapv*dx + self.LagTest*self.Lapu*dx
        self.Lap_Lag_inv = self.Lap_Lag_lhs.Assemble().mat.Inverse(inverse='pardiso')

    def Lap_Inv(self,gfu_rhs,opt):
        if opt == 'peri':
            self.Lap_Lag_rhs = LinearForm(self.fesc_peri_Lag)
            gfu_Lag = GridFunction(self.fesc_peri_Lag)
            difvu, Lag_Mul = gfu_Lag.components
            self.Lap_Lag_rhs += gfu_rhs*self.Lapv*dx
            gfu_Lag.vec.data = self.Lap_Lag_inv*(self.Lap_Lag_rhs.Assemble().vec)
        elif opt == 'Dirichlet':
            self.Lap_rhs = LinearForm(self.fesc_Dirichlet)
            difvu = GridFunction(self.fesc_Dirichlet)
            self.Lap_rhs += gfu_rhs*self.Lapv*dx
            difvu.vec.data = self.Lap_Diri_inv*(self.Lap_rhs.Assemble().vec)
        return difvu

    def GenerateD_mat(self):
        '''
            Gauss points and differential matrix of Lagrangian basis on canonical interval [-1,1]

            nodes: Gauss points in [-1,1]
        '''
        nodes, weights    = GL.generate_Gauss_formula(-1, 1, self.n_collo)
        nodes_with_un     = np.zeros(self.n_collo + 1)
        nodes_with_un[0]  = - 1
        nodes_with_un[1:] = nodes[:]
        D = GL.Diff_Lagrange_Basis_Matrix(nodes_with_un)
        # transformed to y \in [tn,tn+1] : x = 2/tau y - (tn+tn+1)/(tn+1-tn)
        Interval_Change = 2/self.dt
        D = D*Interval_Change
        D_lhs = D[1:,1:]
        D_rhs = D[1:,0]
        # transform numpy matrix to ngsolve
        D_lhs_ng = CF(tuple(D_lhs.flatten()),dims=(self.n_collo,self.n_collo))
        D_rhs_ng = CF(tuple(D_rhs.flatten()),dims=(self.n_collo,1))

        # Compute values of all Lagrange bases at x=1
        extrap_coeff = np.zeros(self.n_collo+1)
        extrap_coeff[0] = GL.Lagrange_Basis(0,nodes_with_un,1)
        for ii in range(self.n_collo):
            extrap_coeff[ii+1] = GL.Lagrange_Basis(ii+1,nodes_with_un,1)
        return D_lhs_ng, D_rhs_ng, extrap_coeff, nodes

    def WeakForm4FixP_Iter(self):
        # After initializing U_N_old
        self.lhs_Niter += 1j*InnerProduct(self.D_mat*self.U_N, self.V_N)*dx - InnerProduct(grad(self.U_N),grad(self.V_N))*dx
        # fold term
        for ii in range(self.n_collo):
            uii = self.U_N_old.components[ii]
            fu2 = self.f_norm2(Norm(uii))
            self.rhs_Niter += -fu2*(uii*self.V_N[ii])*dx
        self.rhs_Niter += -1j*(self.un*InnerProduct(self.D_row,self.V_N))*dx

    def FixPIteration(self):
        '''
            By self.un, set initial data for Newton iteration (U_N_old), each iteration updates U_N_old.

            Obtain self.un_update by extrapolation.
        '''
        err_tmp = GridFunction(self.fesc_working)
        for ii in range(self.n_collo):
            self.U_N_old.components[ii].vec.data = self.un.vec.data
        self.NLS_iters = 0
        err_N = 0
        while err_N>self.NLS_Threshold or self.NLS_iters==0:
            self.lhs_Niter.Assemble()
            self.rhs_Niter.Assemble()
            self.U_N_update.vec.data = self.lhs_Niter.mat.Inverse(self.fes_ALL.FreeDofs(),inverse='pardiso')*self.rhs_Niter.vec
            
            err_N = 0
            for ii in range(self.n_collo):
                err_tmp.vec.data = self.U_N_update.components[ii].vec-self.U_N_old.components[ii].vec
                err_N += Integrate(InnerProduct(err_tmp,err_tmp),self.mesh,element_wise=False)
            err_N = np.sqrt(err_N)
            self.NLS_iters += 1
            self.U_N_old.vec.data = self.U_N_update.vec
        print('Fixed Point Iterations Ends with iters={} and error={}'.format(self.NLS_iters, err_N))

        vec_update = self.extrap_coeff[0]*self.un.vec.FV().NumPy().copy()
        for ii in range(self.n_collo):
            vec_update += self.extrap_coeff[ii+1]*(self.U_N_update.components[ii].vec.FV().NumPy().copy())
        self.un_aux.vec.data = BaseVector(vec_update)
            
    def IniByProj(self):
        pass

    def PP_Pic(self):
        pass

    def GetMassEnergy(self):
        Mass = 1/2*self.L2Omeg(self.un)
        Energy = 1/2*(self.H1semi(self.un) - self.IntOmeg(self.F_norm2(Norm(self.un))))
        return Mass, Energy

    def Projection(self,Eig_opt=False,IO=True):
        '''
            add orthogonal part to ensure the conservation of energy after self.un is updated by iteration method

            *  Eig_opt=True means no projection procedure
        '''
        self.alpha = 1
        self.beta = 0

        if not Eig_opt:
            self.Param_iters = 0
            LapInv = self.Lap_Inv(self.f_norm2(Norm(self.un_aux))*self.un_aux,opt=self.Lap_opt)
            vh = GridFunction(self.fesc_working)
            vh.vec.data = LapInv.vec + self.un_aux.vec

            Lap_inv_u = self.Lap_Inv(self.un_aux,opt=self.Lap_opt)
            Coef = self.IntOmeg(InnerProduct(vh,self.un_aux))/self.IntOmeg(InnerProduct(Lap_inv_u,self.un_aux))
            uh_perp = GridFunction(self.fesc_working)
            uh_perp.vec.data = vh.vec - Coef*Lap_inv_u.vec

            uL2 = self.L2Omeg(self.un_aux)
            uH1_Semi = self.H1semi(self.un_aux)
            uperpL2 = self.L2Omeg(uh_perp)
            uperpH1_Semi = self.H1semi(uh_perp)

            crossprod = self.IntOmeg(InnerProduct(grad(self.un_aux),grad(uh_perp)))
            aux_val = crossprod - self.IntOmeg(self.f_norm2(Norm(self.un_aux))*self.un_aux*Conj(uh_perp))
            gamma = np.arctan2(aux_val.imag,aux_val.real)

            Mass, Energy = self.GetMassEnergy()
            A = 1/2*uL2.real
            B = 1/2*uperpL2.real
            C = 1/2*uH1_Semi.real
            E = 1/2*uperpH1_Semi.real
            D = (np.exp(-1j*gamma)*crossprod).real

            u_interm = GridFunction(self.fesc_working)
            err_alpha_beta = np.inf
            while err_alpha_beta>self.Param_Tolerance and (self.Param_iters<20):
                # each iteration
                u_interm.vec.data = BaseVector(
                    self.alpha*self.un_aux.vec.FV().NumPy()
                    +self.beta*np.exp(1j*gamma)*uh_perp.vec.FV().NumPy()
                )

                Res1 = A*self.alpha**2 + B*self.beta**2 - Mass
                Res2 = C*self.alpha**2 + D*self.alpha*self.beta + E*self.beta**2 - 1/2*self.IntOmeg(self.F_norm2(Norm(u_interm))) - Energy

                Falpha = 1/2*self.IntOmeg(
                    self.f_norm2(Norm(u_interm))*(
                        2*self.alpha*Norm(self.un_aux)**2 +
                        self.beta*(
                            self.un_aux*exp(-1j*gamma)*Conj(uh_perp)+
                            Conj(self.un_aux)*exp(1j*gamma)*uh_perp
                        )
                    )
                )
                Fbeta = 1/2*self.IntOmeg(
                    self.f_norm2(Norm(u_interm))*(
                        2*self.beta*Norm(uh_perp)**2 +
                        self.alpha*(
                            self.un_aux*exp(-1j*gamma)*Conj(uh_perp)+
                            Conj(self.un_aux)*exp(1j*gamma)*uh_perp
                        )
                    )
                )
                F = 1/2*self.IntOmeg(self.F_norm2(Norm(u_interm)))

                Lhs_Mat = np.array([[2*A*self.alpha, 2*B*self.beta], 
                    [2*C*self.alpha+D*self.beta-Falpha, 2*E*self.beta+D*self.alpha-Fbeta]])
                Rhs_Vec = np.array([[Mass + A*self.alpha**2 + B*self.beta**2],
                    [Energy+C*self.alpha**2+D*self.alpha*self.beta+E*self.beta**2+F-Falpha*self.alpha-Fbeta*self.beta]])
                alpha_new, beta_new = np.linalg.solve(Lhs_Mat, Rhs_Vec).flatten()
                err_alpha_beta = np.sqrt(abs(alpha_new-self.alpha)**2 + abs(beta_new-self.beta)**2)
                self.alpha = alpha_new
                self.beta = beta_new
                self.Param_iters += 1
            
            if IO:
                print('Solving Params using iters {}, with alpha, beta = {},{}, err_Alpha_beta = {}'.format(self.Param_iters,self.alpha,self.beta,err_alpha_beta))
            self.un.vec.data = BaseVector(
                self.alpha*self.un_aux.vec.FV().NumPy()
                +self.beta*np.exp(1j*gamma)*uh_perp.vec.FV().NumPy()
            )
        else:
            self.un.vec.data = BaseVector(
                self.alpha*self.un_aux.vec.FV().NumPy()
            )

    def Solving(self):
        while (self.t<self.T) & (abs(self.t-self.T)>1e-10):
            self.t += self.dt
            self.FixPIteration()
            self.Projection()
            self.PP_Pic()

class Coll_Proj_NLS_Newt(Coll_Proj_NLS_FixP):
    def __init__(self, mesh, kappa, p, order, n_collocation, dt, T, thres, quad_order, ref_order) -> None:
        super().__init__(mesh, kappa, p, order, n_collocation, dt, T, thres, quad_order, ref_order)

    def WeakForm4Newton_Iter(self):
        # f expression is kappa x^{(p-1)/2}
        # couple the complex fem space to cover the conjugate part
        self.fes_ALL_db = self.fes_ALL**2
        self.U, self.V = self.fes_ALL_db.TnT()
        self.U_N, self.U_N_conj = self.U[0,:], self.U[1,:]
        self.V_N, self.V_N_conj = self.V[0,:], self.V[1,:]
        self.U_g = GridFunction(self.fes_ALL_db)
        self.U_g_update = GridFunction(self.fes_ALL_db)
        self.U_N_old, self.U_N_old_conj = self.U_g.components[0], self.U_g.components[1]
        self.lhs_Niter = BilinearForm(self.fes_ALL_db)
        self.rhs_Niter = LinearForm(self.fes_ALL_db)
        # After initializing U_N_old
        self.lhs_Niter += 1j*InnerProduct(self.D_mat*self.U_N, self.V_N)*dx 
        self.lhs_Niter += -1j*InnerProduct(self.D_mat*self.U_N_conj, self.V_N_conj)*dx
        self.lhs_Niter += -InnerProduct(grad(self.U),grad(self.V))*dx
        # fold term
        for ii in range(self.n_collo):
            uii = self.U_N_old.components[ii]
            uii_bar = self.U_N_old_conj.components[ii]
            Uii = self.U_N[ii]
            Uii_bar = self.U_N_conj[ii]

            g1u = self.df_norm2(Norm(uii))*(InnerProduct(uii,uii)) + self.f_norm2(Norm(uii))
            g2u = self.df_norm2(Norm(uii))*uii*uii
            g2u_bar = self.df_norm2(Norm(uii))*uii_bar*uii_bar
            fu2 = self.f_norm2(Norm(uii))

            self.rhs_Niter += -fu2*(uii*self.V_N[ii]+uii_bar*self.V_N_conj[ii])*dx
            self.rhs_Niter += (g1u+g2u_bar)*(uii*self.V_N[ii])*dx \
                            + (g2u+g1u)*(uii_bar*self.V_N_conj[ii])*dx
            self.lhs_Niter += (g1u+g2u_bar)*(Uii*self.V_N[ii])*dx \
                            + (g2u+g1u)*(Uii_bar*self.V_N_conj[ii])*dx
        self.rhs_Niter += -1j*(self.un*InnerProduct(self.D_row,self.V_N))*dx\
                            + 1j*(Conj(self.un)*InnerProduct(self.D_row,self.V_N_conj))*dx
        if MatInvOpt == 'Iter':
            self.Pre = Preconditioner(self.lhs_Niter,'multigrid')

    def Newton_Iteration(self,IO=True):
        '''
            By self.un, set initial data for Newton iteration (U_N_old), each iteration updates U_N_old.

            Obtain self.un_update by extrapolation.
        '''
        err_tmp = GridFunction(self.fesc_working)
        for ii in range(self.n_collo):
            self.U_N_old.components[ii].vec.data = self.un.vec
            self.U_N_old_conj.components[ii].vec.data = BaseVector(self.un.vec.FV().NumPy().conj())
        self.NLS_iters = 0
        err_N = np.inf
        # Newton iteration with the given threshold
        # CG_Max_Residue = 0
        # Max_Iters = 0
        N_Iter_CG_Info = GatherData(IO) 
        while err_N>self.NLS_Threshold:
            self.lhs_Niter.Assemble()
            self.rhs_Niter.Assemble()
            if MatInvOpt == 'Direct':
                self.U_g_update.vec.data = self.lhs_Niter.mat.Inverse(self.fes_ALL_db.FreeDofs(),inverse='pardiso')*self.rhs_Niter.vec
            elif MatInvOpt == 'Iter':
                inv = CGSolver(self.lhs_Niter.mat, self.Pre.mat, conjugate=True, maxiter=100, printrates=False,callback=N_Iter_CG_Info.collect)
                self.U_g_update.vec.data = inv*self.rhs_Niter.vec    
            err_N = 0
            # u_aux has many components
            for ii in range(self.n_collo):
                err_tmp.vec.data = self.U_g_update.components[0].components[ii].vec-self.U_N_old.components[ii].vec
                err_N += self.L2Omeg(err_tmp)
            err_N = np.sqrt(err_N)
            self.NLS_iters += 1
            self.U_g.vec.data = self.U_g_update.vec
        if IO:
            N_Iter_CG_Info.print()
            print('Newton Iterations Ends with iters={} and error={}'.format(self.NLS_iters, err_N))

        vec_update = self.extrap_coeff[0]*self.un.vec.FV().NumPy().copy()
        for ii in range(self.n_collo):
            # U_N_old, first component of U_g, has already been updated
            vec_update += self.extrap_coeff[ii+1]*(self.U_N_old.components[ii].vec.FV().NumPy().copy())
        self.un_aux.vec.data = BaseVector(vec_update)

    def Solving(self,save_obj:Vtk_out=None,Eig_opt=False,ThreadNum=10):
        SetNumThreads(ThreadNum)
        with TaskManager():
            while (self.t<self.T) & (abs(self.t-self.T)>1e-10):
                if save_obj is None:
                    pass
                else:
                    save_obj.Output(mesh=self.mesh,tnow=self.t,function=[Norm(self.un)],command='')
                if self.t == 0:
                    self.PP_Pic(True)
                self.t += self.dt
                print('{}: Now time is {}'.format(LogTime(),self.t))
                self.Newton_Iteration()
                self.Projection(Eig_opt)
                self.PP_Pic()
                if save_obj is None:
                    pass
                elif abs(self.t-self.T)<1e-10:
                    save_obj.Output(mesh=self.mesh,tnow=self.t,function=[Norm(self.un)],command='')

class CPN_Converg_1d_Focus_FixP(Coll_Proj_NLS_FixP):
    def __init__(self, mesh, order, n_collocation, dt, T, N_thres) -> None:
        # Model Param
        kappa = 2
        p = 3
        super().__init__(mesh, kappa, p, order, n_collocation, dt, T, N_thres)
        self.tParam = Parameter(0)
        self.exactu_Proj = GridFunction(self.fesc_working)
        self.eu = GridFunction(self.fesc_working)
        self.H1err_set = []
        self.L2err_set = []
        self.mas_c_set = []
        self.eng_c_set = []

    def Sol_Setting(self):
        # related to kappa and p
        self.exact_u_ng = 2/(exp(x+4*self.tParam)+exp(-x-4*self.tParam))*exp(-1j*(2*x+3*self.tParam))

    def IniByProj(self):
        t0 = self.tParam.Get()
        assert t0 == 0
        self.un.Set(self.exact_u_ng)

    def LinfH1Err(self):
        print('L^\infty(0,{},H^1) error is {}'.format(self.T, max(self.H1err_set)))

    def PP_Pic(self):
        # Now: self.un corresponds to un at self.t
        self.tParam.Set(self.t)
        # Projeciton of the exact solution at t in a finite element space (of order ref_order higher)
        self.exactu_Proj.Set(self.exact_u_ng)
        self.eu.vec.data = self.un.vec - self.exactu_Proj.vec
        L2err_instant = np.sqrt(Integrate(InnerProduct(self.eu, self.eu),self.mesh,element_wise=False))
        H1err_instant = np.sqrt(Integrate(InnerProduct(grad(self.eu), grad(self.eu)),self.mesh,element_wise=False)+L2err_instant**2)
        self.L2err_set.append(L2err_instant)
        self.H1err_set.append(H1err_instant)
        Mht, Eht = self.GetMassEnergy()
        self.mas_c_set.append(Mht)
        self.eng_c_set.append(Eht)

class CPN_Converg_1d_Focus_Newt(Coll_Proj_NLS_Newt):
    def __init__(self, mesh, kappa, p, order, n_collocation, dt, T, N_thres, quad_order=5, ref_order=1) -> None:
        # Model Param focusing: kappa = 2, p = 3
        self.kappa = kappa
        self.p = p
        self.PPtime = 0
        super().__init__(mesh, kappa, p, order, n_collocation, dt, T, N_thres, quad_order, ref_order)

    def Sol_Setting(self, exact_u:Callable):
        # Initialize ref sol and numer sol on high order interpolation spaces.
        self.exactu_Proj_ref = GridFunction(self.fesc_working_ref)
        self.u_num_Proj_res = GridFunction(self.fesc_working_ref)
        self.tParam = Parameter(0)
        # Collect err at end points
        self.endt_mas_c_set = []
        self.endt_eng_c_set = []
        self.endt_H1_ex_err_set = []
        self.endt_L2_ex_err_set = []
        self.endt_Newton_it = []
        self.endt_alpbet_it = []
        # Collect err at interior points
        self.int_mas_c_set = []
        self.int_eng_c_set = []
        self.int_H1_ex_err_set = []
        self.int_L2_ex_err_set = []
        self.int_tset = []
        self.int_exactu_H1  = []
        self.endt_exactu_H1 = []
        # Related to kappa and p
        self.exact_u_ng = exact_u(self.tParam)

    def IniByProj(self):
        t0 = self.tParam.Get()
        assert t0 == 0
        self.un.Set(self.exact_u_ng)

    def PostProcess_T(self,t,u_num):
        self.tParam.Set(t)
        # Error by numerical and projection on finer mesh
        self.exactu_Proj_ref.Set(self.exact_u_ng)
        self.u_num_Proj_res.Set(u_num)
        self.eu_ref.vec.data = self.u_num_Proj_res.vec - self.exactu_Proj_ref.vec
        L2_err_instant = np.sqrt(self.L2Omeg(self.eu_ref))
        H1_err_instant = np.sqrt(self.H1Omeg(self.eu_ref))
        Exactu_H1_norm = np.sqrt(self.H1Omeg(self.exactu_Proj_ref))
        # Energy Mass conservation
        Mht = 1/2*self.L2Omeg(u_num)
        Eht = 1/2*self.H1semi(u_num) - 1/2*self.IntOmeg(self.F_norm2(Norm(u_num)))
        return L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm

    def PP_Pic(self,init=False):
        '''
            Caution: self.t should be the end time 
        '''
        if init:
            L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm = self.PostProcess_T(self.t,self.un)
            self.endt_mas_c_set.append(Mht)
            self.endt_eng_c_set.append(Eht)
        else:
            # At this time, self.un corresponds to un, err at end point at self.t
            L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm = self.PostProcess_T(self.t,self.un)
            self.endt_H1_ex_err_set.append(H1_err_instant)
            self.endt_L2_ex_err_set.append(L2_err_instant)
            self.endt_mas_c_set.append(Mht)
            self.endt_eng_c_set.append(Eht)
            self.endt_Newton_it.append(self.NLS_iters)
            self.endt_alpbet_it.append(self.Param_iters)
            self.endt_exactu_H1.append(Exactu_H1_norm)
            # Err at interior collocation points
            for ii,col_x in enumerate(self.collo_inner_nodes):
                t_in_collo = self.t - self.dt + self.dt/2*(col_x+1)
                L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm = self.PostProcess_T(t_in_collo, self.U_N_old.components[ii])
                self.int_H1_ex_err_set.append(H1_err_instant)
                self.int_L2_ex_err_set.append(L2_err_instant)
                self.int_mas_c_set.append(Mht)
                self.int_eng_c_set.append(Eht)
                self.int_tset.append(t_in_collo)
                self.int_exactu_H1.append(Exactu_H1_norm)
        
    def Solving(self,save_obj:Vtk_out=None,Eig_opt=False,IO=False,ThreadNum=10):
        # Updating numerical solution
        myPrinter = Printer(10)
        SetNumThreads(ThreadNum)
        with TaskManager():
            while (self.t<self.T) & (abs(self.t-self.T)>1e-10):
                if save_obj is None:
                    pass
                else:
                    save_obj.Output(mesh=self.mesh,tnow=self.t,function=[Norm(self.un)],command='')
                if self.t == 0:
                    self.PP_Pic(True)
                self.t += self.dt
                myPrinter.PrintDuring(self.t,self.T,T_begin=0) 
                if IO:
                    print('{}: Now time is {}'.format(LogTime(),self.t))
                self.Newton_Iteration(IO)
                start_time = time.time()
                self.Projection(Eig_opt,IO)
                end_time = time.time()
                self.PPtime += end_time - start_time
                self.PP_Pic()
                # Save the final fig
                if save_obj is None:
                    pass
                elif abs(self.t-self.T)<1e-10:
                    save_obj.Output(mesh=self.mesh,tnow=self.t,function=[Norm(self.un)],command='')

class CPN_Converg_2d_Newt(CPN_Converg_1d_Focus_Newt):
    def __init__(self, mesh, kappa, p, order, n_collocation, dt, T, N_thres, ref_order, quad_order=5) -> None:
        super().__init__(mesh, kappa, p, order, n_collocation, dt, T, N_thres, quad_order, ref_order)

class PP_Process:
    def Sol_Setting(self, exact_u:Callable):
        # Initialize ref sol and numer sol on high order interpolation spaces.
        self.exactu_Proj_ref = GridFunction(self.fesc_working_ref)
        self.u_num_Proj_res = GridFunction(self.fesc_working_ref)
        self.tParam = Parameter(0)
        # Collect err at end points
        self.endt_mas_c_set = []
        self.endt_eng_c_set = []
        self.endt_H1_ex_err_set = []
        self.endt_L2_ex_err_set = []
        self.endt_Newton_it = []
        self.endt_alpbet_it = []
        # Collect err at interior points
        self.int_mas_c_set = []
        self.int_eng_c_set = []
        self.int_H1_ex_err_set = []
        self.int_L2_ex_err_set = []
        self.int_tset = []
        self.int_exactu_H1  = []
        self.endt_exactu_H1 = []
        self.exact_u_ng = exact_u(self.tParam)

    def IniByProj(self):
        t0 = self.tParam.Get()
        assert t0 == 0
        self.un.Set(self.exact_u_ng)

    def PostProcess_T(self,t,u_num):
        self.tParam.Set(t)
        # Error by numerical and projection on finer mesh
        self.exactu_Proj_ref.Set(self.exact_u_ng)
        self.u_num_Proj_res.Set(u_num)
        self.eu_ref.vec.data = self.u_num_Proj_res.vec - self.exactu_Proj_ref.vec
        L2_err_instant = np.sqrt(self.L2Omeg(self.eu_ref))
        H1_err_instant = np.sqrt(self.H1Omeg(self.eu_ref))
        Exactu_H1_norm = np.sqrt(self.H1Omeg(self.exactu_Proj_ref))
        # Energy Mass conservation
        Mht = 1/2*self.L2Omeg(u_num)
        Eht = 1/2*self.H1semi(u_num) - 1/2*self.IntOmeg(self.F_norm2(Norm(u_num)))
        return L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm

    def PP_Pic(self,init=False):
        if init:
            L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm = self.PostProcess_T(self.t,self.un)
            self.endt_mas_c_set.append(Mht)
            self.endt_eng_c_set.append(Eht)
        else:
            L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm = self.PostProcess_T(self.t,self.un)
            self.endt_H1_ex_err_set.append(H1_err_instant)
            self.endt_L2_ex_err_set.append(L2_err_instant)
            self.endt_mas_c_set.append(Mht)
            self.endt_eng_c_set.append(Eht)
            self.endt_Newton_it.append(self.NLS_iters)
            self.endt_alpbet_it.append(self.Param_iters)
            self.endt_exactu_H1.append(Exactu_H1_norm)
            # Err at interior collocation points
            for ii,col_x in enumerate(self.collo_inner_nodes):
                t_in_collo = self.t + self.dt/2*(col_x-1)
                L2_err_instant, H1_err_instant, Mht, Eht, Exactu_H1_norm = self.PostProcess_T(t_in_collo, self.U_N_old.components[ii])
                self.int_H1_ex_err_set.append(H1_err_instant)
                self.int_L2_ex_err_set.append(L2_err_instant)
                self.int_mas_c_set.append(Mht)
                self.int_eng_c_set.append(Eht)
                self.int_tset.append(t_in_collo)
                self.int_exactu_H1.append(Exactu_H1_norm)
