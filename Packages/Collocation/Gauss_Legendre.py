import numpy as np
from typing import Union

#========Gauss_weight, Gauss_point==============================================================
def generate_Gauss_formula(vertice_left, vertice_right, Gauss_type):
    '''
        Generate Guass points and weights of Gauss quarature rule on interval [x_l, x_r]. 
        
        For Gauss_type k, algebraic accuracy is 2k-1.
    '''
    h = vertice_right - vertice_left
    ctr = (vertice_right + vertice_left)/2
    
    if Gauss_type == 1:
        Gauss_weight_normal = np.array([2])
        Gauss_point_normal  = np.array([0])

    elif Gauss_type == 2:
        Gauss_weight_normal = np.array([1, 1])
        Gauss_point_normal  = np.array([-0.57735026918962576450914878050195745564760175127012656, 
                                        0.57735026918962576450914878050195745564760175127012656])

    elif Gauss_type == 3:
        Gauss_weight_normal = np.array([0.555555555555555555555555555555555555555555555555555556, 
                                        0.8888888888888888888888888888888888888888888888888888889, 
                                        0.55555555555555555555555555555555555555555555555555556])
        Gauss_point_normal  = np.array([-0.77459666924148337703585307995647992216658434105831767, 0, 
                                        0.77459666924148337703585307995647992216658434105831767])

    elif Gauss_type == 4:
        Gauss_point_normal  = np.array([-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116])
        Gauss_weight_normal = np.array([ 0.3478548461,  0.6521451549, 0.6521451549, 0.3478548461])

    elif Gauss_type == 5:
        Gauss_point_normal  = np.array([-0.9324695142, -0.6612093865, -0.2386191761, 0.2386191761,\
                              0.6612093865,  0.9324695142])
        Gauss_weight_normal = np.array([ 0.1713244924,  0.3607615730,  0.4679139346, 0.4679139346,\
                              0.3607615730,  0.1713244924])
    else: raise Exception('No such Gauss formula')
        
    Gauss_point = Gauss_point_normal*h/2 + ctr
    Gauss_weight = Gauss_weight_normal*h/2
        
    return  Gauss_point, Gauss_weight

#==Generte-Legendre-Polys============================================
def LegendrePoly(n,x:np.ndarray):
    '''
        Compute the n-th Legendre poly's (orthogonal poly on [-1,1]) values (and derivatives) at x.

        nPn(x) = (2n-1)xPn-1(x) - (n-1)Pn-2(x)
    '''
    m = len(x)
    
    polylst = np.ones(m)    # L_0 = 1,
    poly = x                # L_1 = x
    pderlst = np.zeros(m)    # L'_0 = 0
    pder = np.ones(m)        # L'_1 = 1
    
    if n == 0:      
        polyn = polylst
        pdern = pderlst
    elif n == 1:
        polyn = poly
        pdern = pder
    else:
        pdern, polyn = np.zeros(m), np.zeros(m) # L_N vs L'_n 
        for k in range(2,n+1):                             # Three-term recurrence relation:  
            polyn = (2*k-1)*x*poly/k - (k-1)*polylst/k     # kL_k(x)=(2k-1)xL_{k-1}(x)-(k-1)L_{k-2}(x)
            pdern = pderlst+(2*k-1)*poly                   # L_k'(x)=L_{k-2}'(x)+(2k-1)L_{k-1}(x)
            polylst = poly 
            poly = polyn
            pderlst = pder 
            pder = pdern
    return polyn, pdern

#==Generte-Differential-Matrix=======================================
def GaussLegendreDiff(Gauss_type, nodes):
    
    x, n = nodes, Gauss_type
    if Gauss_type == 0: D = []
    D = np.zeros((n, n))    
    y, dy = LegendrePoly(n,x)
    
    for  k in range(0, n):
        for j in range(0, n):
            if k == j: D[k][j] = x[k]/(1 - x[k]**2) 
            else:      D[k][j] = dy[k]/dy[j]/(x[k] - x[j])

    return D
    
#==================================================================
# https://math.stackexchange.com/questions/809927/first-derivative-of-lagrange-polynomial
def Diff_Lag_Basis_on_Nodes(Lag_Basis_j:int, nodes):
    '''
        Value of differential of j-th (in nodes) Lagrangian basis at all nodes 
        
        if node = basis, using Logarithmic derivative

        else: Lebniz rule, most terms are zero.
    '''
    n = len(nodes)
    Diff_L_Basis_on_Nodes = np.zeros(n)
    for ii in range(n):
        if ii == Lag_Basis_j:
            # Logarithmic derivative
            res = 0
            for kk in range(0, n):
                if kk != Lag_Basis_j:
                    res += 1/(nodes[Lag_Basis_j] - nodes[kk])
            Diff_L_Basis_on_Nodes[ii] = res*Lagrange_Basis(Lag_Basis_j, nodes, nodes[ii]) # actually 1
        else:
            res = 1
            for kk in range (0, n):
                if kk != Lag_Basis_j and kk != ii:
                    res = res*(nodes[ii] - nodes[kk])/(nodes[Lag_Basis_j] - nodes[kk])
            Diff_L_Basis_on_Nodes[ii] = res/(nodes[Lag_Basis_j] - nodes[ii])   
    return Diff_L_Basis_on_Nodes

#==================================================================
def Lagrange_Basis(Lag_Basis_j:int,nodes:np.ndarray,x:Union[np.ndarray,float]):
    '''
        Given interpolation nodes, compute the j Lagrangian polynomial at x

        Return ndarray the same size as x
    '''
    n = len(nodes)
    Lj_x = 1
    for k in range (0, n):
        if k != Lag_Basis_j:
            Lj_x = Lj_x*(x - nodes[k])/(nodes[Lag_Basis_j] - nodes[k])
    return Lj_x

#==================================================================
def Diff_Lagrange_Basis_Matrix(nodes):    
    '''
        Lagrangian Basis Differential Matrix --- Interpolate Nodes
        
        Dij: value of j-th basis on i-th node
    '''
    N = len(nodes)
    D = np.zeros((N, N))
    for i in range (0, N):
        # 第i列是第i个基函数在其余点的导数值
        D[:,i] = Diff_Lag_Basis_on_Nodes(i, nodes)
    return D