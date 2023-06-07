from cmath import inf
import os
from netgen.geom2d import *
from ngsolve import *
import numpy as np
from collections import Counter
import scipy
from scipy.sparse import linalg as LA
import triangle as tr
import netgen.meshing as ngm
from pyevtk.hl import polyLinesToVTK
from Package_MyCode import FO
from Package_G1dMesh import Mesh1d

try:
    from Package_MyNgFunc import NGMO
except:
    pass

Rot = lambda x: (x[1],-x[0])
def Rotation2d(x: tuple, theta: float):
    RotMat = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    NewPos = RotMat@np.array(x)
    return (NewPos[0], NewPos[1])

class Vtk_out:
    def __init__(self, T_set_old, n_set, pathname, T_begin=0):
        '''
            Generate to be saved Tsets
        '''
        T_set = [T_i for T_i in T_set_old if T_i>T_begin]
        n_set = [n_set[ii] for ii in range(len(n_set)) if T_set_old[ii]>T_begin]
        if len(T_set)==1:
            self.Tsets = np.linspace(T_begin,T_set[0],n_set[0]+1)
        else:
            self.Tsets = np.zeros(1)
            for ii in range(len(T_set)):
                if ii == 0:
                    init,endt = T_begin,T_set[0]
                else:
                    init,endt = T_set[ii-1],T_set[ii]
                self.Tsets = np.append(self.Tsets,np.linspace(init,endt,n_set[ii]+1)[1:])

        self.index = 0
        self.fillnum = len(str(max(n_set)))
        self.pathname = pathname
        self.done = False
        if not os.path.exists(pathname):
            os.makedirs(pathname)
            mystr = '  vtk files are saved in {}  '.format(pathname)
            print('{:#^60}'.format(mystr))
        else:
            self.done = True
        self.Rel_Mapping = {}
        np.save(file=os.path.join(pathname,'Tset.npy'),arr=self.Tsets)
    
    def GenerateMapping(self, filepath, LoadT=0):
        '''
            LoadT = 0 version
        '''
        vtu_list = [fname for fname in FO.Get_File_List(filepath) if fname.endswith('vtu')]
        n_file = len(vtu_list)
        myTSets = [t for t in self.Tsets if t>=LoadT]
        for ii, T in enumerate(myTSets):
            if ii<n_file:
                self.Rel_Mapping[vtu_list[ii].split('.vtu')[0]] = T
        np.save(file=os.path.join(filepath,'Rel_Mapping.npy'),arr=self.Rel_Mapping,allow_pickle=True)

    def Output(self, mesh, function, tnow, command='', names=['sol'],subdivision=0):
        perform = False
        if command == 'do':
            perform = True
        elif tnow >= self.Tsets[self.index]:
            perform = True
        if perform:
            vtk = VTKOutput(ma=mesh,coefs=function,names=names,\
                            filename=os.path.join(self.pathname,('{}'.format(self.index)).zfill(self.fillnum)),\
                            subdivision=subdivision,legacy=False)
            vtk.Do()
            self.index += 1

    def Generate_PVD(self,pvd_name):
        # generate only one path
        FO.PVD_Generate(pvd_path=self.pathname,folder_path_set=[''], pvd_name=pvd_name)
class Vtk_out_BND(Vtk_out):
    def __init__(self, T_set, n_set, pathname, T_begin=0):
        super().__init__(T_set, n_set, pathname, T_begin)
    
    def Output(self, mesh, function, tnow, command='', names=['sol']):
        perform = False
        if command == 'do':
            perform = True
        elif tnow >= self.Tsets[self.index]:
            perform = True
        if perform:
            local_name = ('{}'.format(self.index)).zfill(self.fillnum)
            file_path = os.path.join(self.pathname,local_name)
            vtk = VTKOutput(ma=mesh,coefs=function,names=names,\
                            filename=file_path,\
                            subdivision=0,legacy=False)
            self.Rel_Mapping[local_name] = tnow
            vtk.Do(vb=BND)
            self.index += 1
        return perform

class Vtk_out_1d(Vtk_out):
    '''
        Save vtk file of 1d curve in plane (z=0) by coordinates

        Example: 
            # generate vtk output object
            T     = [1]
            n_vtk = [10]
            VTU_Path = './vtkfile'
            Coords2d0 = np.array([[np.cos(theta), np.sin(theta)] 
                        for theta in np.linspace(0,2*np.pi,20)[:-1]])
            vtk_Obj_1d = Vtk_out_1d(T,n_vtk,VTU_Path)
            for t in np.linspace(0,1,20):
                # radius changes from 2 to 1
                Coords2d = (2-t)*Coords2d0
                perform_res = vtk_Obj_1d.Output(Coords2d,None,tnow=t)
            np.save(file=os.path.join(VTU_Path,'Rel_Mapping.npy'),arr=vtk_Obj_1d.Rel_Mapping,allow_pickle=True)
            vtk_Obj_1d.Generate_PVD('test.pvd')
    '''
    
    def __init__(self, T, n, pathname, T_begin=0):
        super().__init__(T, n, pathname, T_begin)
        self.VertsCoords = {}
    
    def Output(self, vertices, pData_dict=None, tnow=None):
        '''
            输入一维曲线的逆时针序节点，输出vtk文件
        '''
        perform = False
        if pData_dict is None:
            pData_dict = {'zero': np.zeros(vertices.shape[0])}
        for key,val in pData_dict.items():
            pData_dict[key] = np.append(val,val[0])
        if self.index<len(self.Tsets) and (tnow is None or tnow >= self.Tsets[self.index]):
            # Positions of points that define lines
            npoints = len(vertices)
            x = np.zeros(npoints+1)
            y = np.zeros(npoints+1)
            z = np.zeros(npoints+1)

            # First line
            for ii, v in enumerate(vertices):
                x[ii], y[ii], z[ii] = vertices[ii,0], vertices[ii,1], 0.0
            x[-1], y[-1], z[-1] = vertices[0,0], vertices[0,1], 0.0

            # Connectivity of the lines
            pointsPerLine = np.zeros(1)
            pointsPerLine[0] = len(vertices)+1

            local_name = ('{}'.format(self.index)).zfill(self.fillnum)
            file_path = os.path.join(self.pathname,local_name)
            polyLinesToVTK(
                file_path,
                x,
                y,
                z,
                pointsPerLine=pointsPerLine,
                pointData=pData_dict
            )
            self.Rel_Mapping[local_name] = tnow
            self.index += 1
            perform = True
        return perform
        
    def LineSave(self, vertices, tnow=None):
        if self.index<len(self.Tsets) and (tnow is None or tnow >= self.Tsets[self.index]):
            self.VertsCoords[tnow] = (vertices)
            self.index += 1

    def SaveNpy(self,pvd_name):
        np.save(os.path.join(self.pathname,pvd_name), self.VertsCoords, allow_pickle=True)

class yxt_1d_out(Vtk_out):
    def __init__(self, x_save_ref, T, n, pathname, T_begin=0):
        super().__init__(T, n, pathname, T_begin)
        self.Data_List = []
        self.t_List = []
        self.x_ref = x_save_ref
        Mesh_Obj = Mesh1d(min(x_save_ref),max(x_save_ref),len(x_save_ref)-1)
        self.mesh     = Mesh(Mesh_Obj.ngmesh)
        self.save_fem = H1(self.mesh)
        self.funvalue = GridFunction(self.save_fem)

    def Output(self, mesh, function, tnow, command='', names=['sol'], subdivision=0):
        assert(len(function)==1)
        func = function[0]
        perform = False
        if command == 'do':
            perform = True
        elif tnow >= self.Tsets[self.index]:
            perform = True
        if perform:
            self.funvalue.Set(func)
            self.Data_List.append(self.funvalue.vec.FV().NumPy().copy())
            self.t_List.append(tnow)
            self.index += 1
        
        

class SplineSeg3():
    '''
        输入三个控制点的NURBF三次样条，特定的weight使得它可以表示ellipic曲线
    '''
    
    def __init__(self,p1,p2,p3) -> None:
        self.p1, self.p2, self.p3 = p1, p2, p3
        d13 = np.sqrt(np.sum((self.p1-self.p3)**2))
        d12 = np.sqrt(np.sum((self.p1-self.p2)**2))
        d23 = np.sqrt(np.sum((self.p2-self.p3)**2))
        self.weight = np.sqrt((2*d13**2)/(d12**2+d23**2))

    def GetPoint(self,t):
        b1, b2, b3 = (1-t)**2, self.weight*t*(1-t), t**2
        hp = b1*self.p1 + b2*self.p2 + b3*self.p3
        w = b1 + b2 + b3
        return (1/w)*hp

    def GetDerivatives(self,t):
        b1 = (1-t)**2
        b2 = self.weight*t*(1-t)
        b3 = t**2
        w = b1+b2+b3
        b1 *= 1/w
        b2 *= 1/w
        b3 *= 1/w

        b1p = 2*(t-1)
        b2p = self.weight*(1-2*t)
        b3p = 2*t
        wp = b1p + b2p + b3p
        fac1 = wp/w
        b1p *= 1/w
        b2p *= 1/w
        b3p *= 1/w

        b1pp = 2
        b2pp = -2*self.weight
        b3pp = 2
        wpp = b1pp + b2pp + b3pp
        fac2 = (wpp*w - 2*wp*wp)/w**2

        point = b1*self.p1 + b2*self.p2 + b3*self.p3

        first = (b1p-b1*fac1) * self.p1 + (b2p - b2*fac1) * self.p2 + (b3p - b3*fac1) * self.p3
        second = (b1pp/w - 2*b1p*fac1 - b1*fac2) * self.p1 + \
            (b2pp/w - 2*b2p*fac1 - b2*fac2) * self.p2 + \
            (b3pp/w - 2*b3p*fac1 - b3*fac2) * self.p3
        return point, first, second
    
    def Project(self,point,break_tag=False):
        t_old = -1
        tmin = 1
        dist_min2 = np.linalg.norm(self.GetPoint(tmin) - point)
        
        # 初始tmin是1，进一步和0，1上四等分的时间点的距离进行比较
        for ti in np.linspace(0,0.99,4):
            di = np.linalg.norm(self.GetPoint(ti) - point)
            if di<dist_min2:
                tmin = ti
                dist_min2 = di
        t = tmin
        i = 0
        # 终止条件分别是：超出了三次样条的参数范围，迭代次数过多以及收敛
        while (t>-0.5 and t<1.5 and i<20 and abs(t-t_old) > 1e-8):
            phi, phip, phipp = self.GetDerivatives(t)
            t_old = t
            phimp = phi-point
            t -= (np.inner(phip,phimp))/(np.inner(phipp,phimp) + np.inner(phip,phip))
            i += 1

        # if (i<20 and t > -0.4 and t < 1.4):
        if True:
            if t<0:
                t = 0
            if t>1:
                t = 1
            
            point_on_curve = self.GetPoint(t)
            dist = np.linalg.norm(point-point_on_curve)
            phi = self.GetPoint(0)
            auxdist = np.linalg.norm(phi-point)
            # 起点的距离更近
            if auxdist<dist:
                t = 0
                point_on_curve = phi
                dist = auxdist
            phi = self.GetPoint(1)
            auxdist = np.linalg.norm(phi-point)
            if auxdist<dist:
                t = 1
                point_on_curve = phi
                dist = auxdist
        else:
            t0 = 0
            t1 = 0.5
            t2 = 1
            while (t2-t0>1e-8):
                phi = self.GetPoint(t0)
                d0 = np.linalg.norm(phi - point)
                phi = self.GetPoint(t1)
                d1 = np.linalg.norm(phi - point)
                phi = self.GetPoint(t2)
                d2 = np.linalg.norm(phi - point)

                a = (2*d0 - 4*d1 + 2*d2)/np.power(t2-t0,2)
                if a<=0:
                    if d0<d2:
                        t2 -= 0.3*(t2-t0)
                    else:
                        t0 += 0.3*(t2-t0)
                    t1 = 0.5*(t2+t0)
                else:
                    b = (d1-d0-a*(t1*t1-t0*t0))/(t1-t0)
                    auxt1 = -0.5*b/a
                    if auxt1<t0:
                        t2 -= 0.4*(t2-t0)
                        t0 = max(0,t0-0.1*(t2-t0))
                    elif auxt1 > t2:
                        t0 += 0.4*(t2-t0)
                        t2 = min(1,t2+0.1*(t2-t0))
                    else:
                        t1 = auxt1
                        auxt1 = 0.25*(t2-t0)
                        t0 = max(0,t1-auxt1)
                        t2 = min(1,t1+auxt1)
                    t1 = 0.5*(t2+t0)

            phi = self.GetPoint(t0)
            d0 = np.linalg.norm(phi - point)   
            phi = self.GetPoint(t1)
            d1 = np.linalg.norm(phi - point)
            phi = self.GetPoint(t2)
            d2 = np.linalg.norm(phi - point)
            mind = d0
            t = t0
            if d1<mind:
                t = t1
                mind = d1
            if d2<mind:
                t = t2
                mind = d2
            point_on_curve = self.GetPoint(t)
        
        return point_on_curve

def PieceSplineSegs(ctr_point):
    '''
        control points of piecewise cubic splines, nx2 ndarray
        Output: list that constructs SplineSeg3
    '''
    N = len(ctr_point)
    ctr_point = np.array(ctr_point)
    ctr_point_set = []
    ii = 0
    while ii<N:
        p0 = ctr_point[ii]
        p1 = ctr_point[ii+1]
        if ii+2 == N:
            p2 = ctr_point[0]
        else:
            p2 = ctr_point[ii+2]
        ctr_point_set.append((p0,p1,p2))
        ii += 2
    Geometry_in_Splines = [SplineSeg3(tripoints[0],tripoints[1],tripoints[2]) for tripoints in ctr_point_set]
    return Geometry_in_Splines

def EllipseInRect(msize,bndname='bnd'):
    geo = SplineGeometry()
    geo.AddRectangle((-1,-1), (1,1),
                     bcs=["b","r","t","l"],
                     leftdomain=1)
    ctr_point = [(2,0), (2,1), (0,1), (-2,1), (-2,0), 
                (-2,-1), (0,-1), (2,-1)]
    ctr_point = [(0.25*x, 0.25*y) for x,y in ctr_point]
    p1,p2,p3,p4,p5,p6,p7,p8 = [ geo.AppendPoint(x,y) for x,y in ctr_point]
    geo.Append (["spline3", p1, p2, p3],bc=bndname,
                     leftdomain=2, rightdomain=1)
    geo.Append (["spline3", p3, p4, p5],bc=bndname,
                     leftdomain=2, rightdomain=1)
    geo.Append (["spline3", p5, p6, p7],bc=bndname,
                     leftdomain=2, rightdomain=1)
    geo.Append (["spline3", p7, p8, p1],bc=bndname,
                     leftdomain=2, rightdomain=1)
    geo.SetMaterial (1, "outer")
    geo.SetMaterial (2, "inner")
    mesh = Mesh(geo.GenerateMesh(maxh=msize))
    return mesh, ctr_point

def CircleInRect(msize,bndname='bnd',radi=0.5, hInterface=None):
    if hInterface is not None:
        pass
    else:
        hInterface = msize
    geo = SplineGeometry()
    geo.AddRectangle((-1,-1), (1,1),
                     bcs=["b","r","t","l"],
                     leftdomain=1)
    ctr_point = [(1,0), (1,1), (0,1), (-1,1), (-1,0), 
                (-1,-1), (0,-1), (1,-1)]
    ctr_point = [(radi*x, radi*y) for x,y in ctr_point]
    p1,p2,p3,p4,p5,p6,p7,p8 = [ geo.AppendPoint(x,y) for x,y in ctr_point]
    geo.Append (["spline3", p1, p2, p3],bc=bndname,
                     leftdomain=2, rightdomain=1, maxh=hInterface)
    geo.Append (["spline3", p3, p4, p5],bc=bndname,
                     leftdomain=2, rightdomain=1, maxh=hInterface)
    geo.Append (["spline3", p5, p6, p7],bc=bndname,
                     leftdomain=2, rightdomain=1, maxh=hInterface)
    geo.Append (["spline3", p7, p8, p1],bc=bndname,
                     leftdomain=2, rightdomain=1, maxh=hInterface)
    geo.SetMaterial (1, "outer")
    geo.SetMaterial (2, "inner")
    mesh = Mesh(geo.GenerateMesh(maxh=msize))
    return mesh, ctr_point

def MoreEllipseInRect(msize,bndname='bnd'):
    geo = SplineGeometry()
    geo.AddRectangle((-1,-1), (1,1),
                     bcs=["b","r","t","l"],
                     leftdomain=1)
    ellipse_set = [{"center":(-0.4,0.4), "axis":(0.3,0.2), "domainIndex":(2,1), "name": "left"},
                   {"center":(0.4,-0.4), "axis":(0.2,0.3), "domainIndex":(2,1), "name": "right"}]
    ctr_point_set = {}
    for ellipse_param in ellipse_set:
        l_axis, s_axis = ellipse_param["axis"]
        l_index, r_index = ellipse_param["domainIndex"]
        ctr_point = [(l_axis,0), (l_axis,s_axis), (0,s_axis), (-l_axis,s_axis), (-l_axis,0), 
                    (-l_axis,-s_axis), (0,-s_axis), (l_axis,-s_axis)]
        cx, cy = ellipse_param["center"]
        ctr_point = [(cx+pos[0], cy+pos[1]) for pos in ctr_point]
        p1,p2,p3,p4,p5,p6,p7,p8 = [geo.AppendPoint(x,y) for x,y in ctr_point]
        geo.Append (["spline3", p1, p2, p3],bc=bndname,
                        leftdomain=l_index, rightdomain=r_index)
        geo.Append (["spline3", p3, p4, p5],bc=bndname,
                        leftdomain=l_index, rightdomain=r_index)
        geo.Append (["spline3", p5, p6, p7],bc=bndname,
                        leftdomain=l_index, rightdomain=r_index)
        geo.Append (["spline3", p7, p8, p1],bc=bndname,
                        leftdomain=l_index, rightdomain=r_index)
        ctr_point_set[ellipse_param["name"]] = ctr_point
    mesh = Mesh(geo.GenerateMesh(maxh=msize))
    return mesh, ctr_point_set

def Ellipse(msize,bndname='bnd'):
    geo = SplineGeometry()
    ctr_point = [(2,0), (2,1), (0,1), (-2,1), (-2,0), 
                (-2,-1), (0,-1), (2,-1)]
    ctr_point = [(0.25*x, 0.25*y) for x,y in ctr_point]
    p1,p2,p3,p4,p5,p6,p7,p8 = [ geo.AppendPoint(x,y) for x,y in ctr_point]
    geo.Append (["spline3", p1, p2, p3],bc=bndname,
                     leftdomain=1, rightdomain=0)
    geo.Append (["spline3", p3, p4, p5],bc=bndname,
                     leftdomain=1, rightdomain=0)
    geo.Append (["spline3", p5, p6, p7],bc=bndname,
                     leftdomain=1, rightdomain=0)
    geo.Append (["spline3", p7, p8, p1],bc=bndname,
                     leftdomain=1, rightdomain=0)
    # geo.SetMaterial (1, "outer")
    mesh = Mesh(geo.GenerateMesh(maxh=msize))
    return mesh, ctr_point

def Circle(msize,bndname='bnd'):
    geo = SplineGeometry()
    geo.AddCircle(c=(0,0),
                r=1,
                bc=bndname,
                leftdomain=1,
                rightdomain=0)
    geo.SetMaterial(1, "inner")
    mesh = Mesh(geo.GenerateMesh(maxh=msize))
    return mesh

def QuatrefoilInRect(msize,bndname='bnd',hInterface=None,csafty=1,rin=0.1,theta=0):
    if hInterface is not None:
        pass
    else:
        hInterface = msize
    ctr_point = [(0.,  0.6), (rin, 0.6), (rin, 0.4), 
            (rin, rin), (0.4, rin), (0.6, rin)]
    ctr_point = [Rotation2d(v, theta) for v in ctr_point]
    ctr_point += [Rot(v) for v in ctr_point]
    ctr_point += [Rot(Rot(v)) for v in ctr_point]

    geo = SplineGeometry()
    plist = [geo.AppendPoint(*pnt) for pnt in ctr_point]
    curves = [[["spline3",plist[i],plist[i+1],plist[i+2]], bndname] for i in range(len(plist)-2) if i%2==0]
    curves.append([["spline3",plist[-2],plist[-1],plist[0]], bndname])
    geo.AddRectangle((-1,-1), (1,1),
                     bcs=["b","r","t","l"],
                     leftdomain=1,rightdomain=0)
    [geo.Append(c,bc=bc,leftdomain=1,rightdomain=2,maxh=hInterface) for c,bc in curves]
    geo.SetMaterial (1, "outer")
    geo.SetMaterial (2, "inner")
    mesh = Mesh(geo.GenerateMesh(maxh=msize,curvaturesafety=csafty))
    return mesh, ctr_point

def Quatrefoil(msize,bndname='bnd',hInterface=None,rin=0.1):
    if hInterface is not None:
        pass
    else:
        hInterface = msize
    Rot = lambda x: (x[1],-x[0])
    ## clockwise points
    ctr_point = [(0.,  0.6), (rin, 0.6), (rin, 0.4), 
            (rin, rin), (0.4, rin), (0.6, rin)]
    ctr_point += [Rot(v) for v in ctr_point]
    ctr_point += [Rot(Rot(v)) for v in ctr_point]

    geo = SplineGeometry()
    plist = [geo.AppendPoint(*pnt) for pnt in ctr_point]
    curves = [[["spline3",plist[i],plist[i+1],plist[i+2]], bndname] for i in range(len(plist)-2) if i%2==0]
    curves.append([["spline3",plist[-2],plist[-1],plist[0]], bndname])
    [geo.Append(c,bc=bc,leftdomain=0,rightdomain=1,maxh=hInterface) for c,bc in curves]
    mesh = Mesh(geo.GenerateMesh(maxh=msize))
    return mesh, ctr_point

def GetInterfaceCoord(fes,mesh,InterfaceIndex='bnd'):
    # Return coordinates at the boundary of mesh indexed by InterfaceIndex
    E_ind = []
    for el in fes.Elements(BND):
        if el.mat in InterfaceIndex:
            for v in el.vertices:
                if v not in E_ind:
                    E_ind.append(v) 
    E_coord = []
    for v in mesh.vertices:
        if v in E_ind:
            E_coord.append(v.point)
    E_coord = np.array(E_coord)
    return E_coord, E_ind
    
def Rot2d(theta):
    '''
        Rotated theta anticlockwise    
    '''
    return np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])

def Distant(Geometry_in_Splines,Coords,pnums):
    '''
        Generate sampling points on the boundary of number pnums
        and compute the nearest distance of all Coords
        input: list of Geometry_in_Splines: SplineSeg3()
    '''
    Ref_Points = []
    t_sets = np.linspace(0,1,pnums)
    for Spline3 in Geometry_in_Splines:
        for t in t_sets:
            Ref_Points.append(Spline3.GetPoint(t))
    Ref_Points = np.array(Ref_Points)
    dist_set = []
    for coord in Coords:
        # 对每一个待求距离的坐标点
        Dists = np.linalg.norm(Ref_Points - coord[None,:],axis=1)
        dist = np.min(Dists)
        dist_set.append(dist)
    return dist_set

def Pullback(omega, t, point_set: np.ndarray, Geometry_in_Splines):
    '''
        将point_set中的点正交投影到Geometry_in_Splines逆时针旋转 omega*t 的位置上
        omega: 逆时针旋转角速度
        t: 旋转时间
        point_set: 待投影的点
        Profile: t=0时的参数化曲线
    '''
    phi = omega*t
    Rot_point_set = (Rot2d(-phi)@point_set.transpose()).transpose()

    point_on_spline_set = np.zeros(Rot_point_set.shape)
    ii = 0
    for point in Rot_point_set:
        mind = np.inf
        minp = np.zeros(2)
        for spline in Geometry_in_Splines:
            point_on_spline = spline.Project(point)
            dist = np.linalg.norm(point_on_spline - point)
            if dist<mind:
                mind = dist
                minp = point_on_spline
        point_on_spline_set[ii] = minp
        ii += 1

    point_on_spline_set = (Rot2d(phi)@point_on_spline_set.transpose()).transpose()
    return point_on_spline_set

def Gen_Geometry_in_Splines3(ctr_point):
    N = len(ctr_point)
    ctr_point = np.array(ctr_point)
    ctr_point_set = []
    ii = 0
    while ii<N:
        p0 = ctr_point[ii]
        p1 = ctr_point[ii+1]
        if ii+2 == N:
            p2 = ctr_point[0]
        else:
            p2 = ctr_point[ii+2]
        ctr_point_set.append((p0,p1,p2))
        ii += 2

    Geometry_in_Splines = [SplineSeg3(tripoints[0],tripoints[1],tripoints[2]) 
                        for tripoints in ctr_point_set]
    return Geometry_in_Splines

class Extension():
    def __init__(self,mesh,Ename) -> None:
        self.Type = Ename
        self.mesh = mesh
        self.fesall = VectorH1(self.mesh,order=1,dirichlet="bnd|b|r|t|l")
        self.Disp = GridFunction(self.fesall)
        # 在位置延拓的情形下，Disp需要用当前位置减去初始位置，因此引入初始位置Initial_Pos
        self.Initial_Pos = GridFunction(self.fesall)
        self.Initial_Pos.Interpolate(CoefficientFunction((x,y)))
        # Biharmonic求解时空间是高次，之后在下面的标量线性空间上投影
        self.fesscalar = H1(self.mesh,order=1)
        self.fesvector = VectorH1(self.mesh,order=1)
        self.ScalarTmp = GridFunction(self.fesscalar)
        self.VectorTmp = GridFunction(self.fesvector)

    def FEMSetting(self):
        pass

    def WeakForm(self):
        pass

    def SetEname(self, myname):
        if myname == 'X':
            self.Ename = 'X'
            print('Position Extension!!')
        elif myname == 'V':
            self.Ename = 'V'
            print('Velocity Extension!!')

class HarmonicExt(Extension):
    def __init__(self, mesh, Ename='X') -> None:
        super().__init__(mesh, Ename)
        print('Harmonic Extension is used!!')

    def FEMSetting(self):
        # FEM需要的weak form, Dirichlet bnd
        self.sol = GridFunction(self.fesall)
        self.gInterface = GridFunction(self.fesall)
        self.gFixBnd = GridFunction(self.fesall)
        self.g0 = GridFunction(self.fesall)
        self.lhs, self.rhs = self.WeakForm()

    def WeakForm(self):
        # 延拓position
        U,V = self.fesall.TnT()
        lhs = BilinearForm(self.fesall)
        lhs += -InnerProduct(grad(U),grad(V))*dx
        rhs = LinearForm(self.fesall)
        # g0 stands for D氏边界条件
        rhs += InnerProduct(grad(self.g0),grad(V))*dx
        return lhs, rhs

    def Solving(self,Pos_interface,dt):
        ## Set Dirichlet boundary condition: By values on given Interface to solve Harmonic Extension
        self.gInterface.Interpolate(Pos_interface,definedon=self.mesh.Boundaries('bnd'))
        if self.Ename == 'X':
            self.gFixBnd.Interpolate(CoefficientFunction((x,y)),definedon=self.mesh.Boundaries('b|r|t|l'))
        elif self.Ename == 'V':
            self.gFixBnd.Interpolate(CoefficientFunction((0,0)),definedon=self.mesh.Boundaries('b|r|t|l'))
        self.g0.vec.data = self.gInterface.vec + self.gFixBnd.vec

        ## solving extension equation
        self.lhs.Assemble()
        self.rhs.Assemble()
        self.sol.vec.data = self.g0.vec+self.lhs.mat.Inverse(freedofs=self.fesall.FreeDofs())*self.rhs.vec
        if self.Ename == 'X':
            ## at this time, sol is the extended position
            self.Disp.vec.data = self.sol.vec - self.Initial_Pos.vec 
        elif self.Ename == 'V':
            ## at this time, sol is the extended velocity
            self.Disp.vec.data += BaseVector(dt*self.sol.vec.FV().NumPy())
            
class V_BiHarmonicExt(Extension):
    '''
        Velocity Biharmonic Extension
    '''
    def __init__(self, mesh, Ename='V') -> None:
        super().__init__(mesh, Ename)
        print('BiHarmonic Extension is used!!')
        
    def FEMSetting(self):
        ## Scalar Biharmonic Extension By Mixed FEM
        # FEM weak form, Dirichlet bnd
        self.n = specialcf.normal(2)
        order = 2
        self.V = HDivDiv(self.mesh, order=order-1)
        self.Q = H1(self.mesh, order=order, dirichlet="bnd|b|r|t|l")
        self.X = FESpace([self.V,self.Q])
        self.sol = GridFunction(self.X)
        self.gInterface = GridFunction(self.X)
        self.gFixBnd = GridFunction(self.X)
        self.g0 = GridFunction(self.X)
        self.scalar_lhs, self.scalar_rhs = self.WeakForm()

    def WeakForm(self):
        sigma, w = self.X.TrialFunction()
        tau, v = self.X.TestFunction()
        tang = lambda u: u-(u*self.n)*self.n
        lhs = BilinearForm(self.X, symmetric=True)
        lhs += (InnerProduct (sigma, tau) + div(sigma)*grad(v) \
            + div(tau)*grad(w) - 1e-10*w*v )*dx \
            + (-(sigma*self.n) * tang(grad(v)) - (tau*self.n)*tang(grad(w)))*dx(element_boundary=True)
        rhs = LinearForm(self.X)
        rhs += (-div(tau)*grad(self.g0.components[1]) + 1e-10*self.g0.components[1]*v)*dx
        rhs += (tau*self.n)*tang(grad(self.g0.components[1]))*dx(element_boundary=True)
        return lhs, rhs
    
    def Solving(self,V_interface,dt):
        ## Set Dirichlet boundary condition: Set values on Interface to solve BiHarmonic Extension
        ## split as two scalar equations
        ## Thus, boundary value g0 is a funciton on [X,Q], thus need to revise [1]component as velocity extension
        self.gInterface.components[1].Interpolate(V_interface[0],definedon=self.mesh.Boundaries('bnd'))
        self.gFixBnd.components[1].Interpolate(CoefficientFunction(0),definedon=self.mesh.Boundaries('b|r|t|l'))
        self.g0.vec.data = self.gInterface.vec + self.gFixBnd.vec
        
        self.scalar_lhs.Assemble()
        self.scalar_rhs.Assemble()

        self.sol.vec.data = self.scalar_lhs.mat.Inverse(self.X.FreeDofs()) * self.scalar_rhs.vec + self.g0.vec
        self.ScalarTmp.Interpolate(self.sol.components[1])
        self.Disp.components[0].vec.data += BaseVector(dt*self.ScalarTmp.vec.FV().NumPy())
        
        ## extend y of velocity
        self.gInterface.components[1].Interpolate(V_interface[1],definedon=self.mesh.Boundaries('bnd'))
        self.gFixBnd.components[1].Interpolate(CoefficientFunction(0),definedon=self.mesh.Boundaries('b|r|t|l'))
        self.g0.vec.data = self.gInterface.vec + self.gFixBnd.vec

        self.scalar_lhs.Assemble()
        self.scalar_rhs.Assemble()

        self.sol.vec.data = self.scalar_lhs.mat.Inverse(self.X.FreeDofs()) * self.scalar_rhs.vec + self.g0.vec
        self.ScalarTmp.Interpolate(self.sol.components[1])
        self.Disp.components[1].vec.data += BaseVector(dt*self.ScalarTmp.vec.FV().NumPy())

class V_BiHarmonicExt2(Extension):
    '''
        Introduce w = Lap u to solve the BiHarmonic equation
    '''
    def __init__(self, mesh, Ename='V') -> None:
        super().__init__(mesh, Ename)
        print('BiHarmonic Extension Method2 is used!!')
        
    def FEMSetting(self):
        ## Scalar Biharmonic Extension By Mixed FEM
        # (w,wt)+(grad u, grad wt) = 0, wt in H1
        # (grad w, grad ut)        = 0, ut in H01
        order = 2
        self.V = H1(self.mesh, order=order)
        self.Q = H1(self.mesh, order=order, dirichlet="bnd|b|r|t|l")
        self.X = FESpace([self.V,self.Q])
        self.sol = GridFunction(self.X)
        self.gInterface = GridFunction(self.X)
        self.gFixBnd = GridFunction(self.X)
        self.g0 = GridFunction(self.X)
        self.scalar_lhs, self.scalar_rhs = self.WeakForm()

    def WeakForm(self):
        # (w,wt)+(grad u, grad wt) = 0, wt in H1
        # (grad w, grad ut)        = 0, ut in H01
        w, u = self.X.TrialFunction()
        wt, ut = self.X.TestFunction()
        lhs = BilinearForm(self.X, symmetric=True)
        lhs += (w*wt + InnerProduct(grad(u),grad(wt)))*dx
        lhs += InnerProduct(grad(w),grad(ut))*dx
        rhs = LinearForm(self.X)
        rhs += -InnerProduct(grad(self.g0.components[1]),grad(wt))*dx
        return lhs, rhs
    
    def Solving(self,V_interface,dt):
        self.gInterface.components[1].Interpolate(V_interface[0],definedon=self.mesh.Boundaries('bnd'))
        self.gFixBnd.components[1].Interpolate(CoefficientFunction(0),definedon=self.mesh.Boundaries('b|r|t|l'))
        self.g0.vec.data = self.gInterface.vec + self.gFixBnd.vec
        
        self.scalar_lhs.Assemble()
        self.scalar_rhs.Assemble()

        self.sol.vec.data = self.scalar_lhs.mat.Inverse(self.X.FreeDofs()) * self.scalar_rhs.vec + self.g0.vec
        self.ScalarTmp.Interpolate(self.sol.components[1])
        self.Disp.components[0].vec.data += BaseVector(dt*self.ScalarTmp.vec.FV().NumPy())
        
        self.gInterface.components[1].Interpolate(V_interface[1],definedon=self.mesh.Boundaries('bnd'))
        self.gFixBnd.components[1].Interpolate(CoefficientFunction(0),definedon=self.mesh.Boundaries('b|r|t|l'))
        self.g0.vec.data = self.gInterface.vec + self.gFixBnd.vec

        self.scalar_lhs.Assemble()
        self.scalar_rhs.Assemble()

        self.sol.vec.data = self.scalar_lhs.mat.Inverse(self.X.FreeDofs()) * self.scalar_rhs.vec + self.g0.vec
        self.ScalarTmp.Interpolate(self.sol.components[1])
        self.Disp.components[1].vec.data += BaseVector(dt*self.ScalarTmp.vec.FV().NumPy())

class X_BiHarmonicExt2(Extension):
    '''
        使用位置来进行双调和延拓： Lap^2 X = 0，边界条件为Dirichlet边界条件，以及 p X/p n = n
    '''
    def __init__(self, mesh, Ename='V') -> None:
        super().__init__(mesh, Ename)
        print('BiHarmonic Extension Method2 is used!!')

    def FEMSetting(self):
        ## Scalar Biharmonic Extension By Mixed FEM, (vector valued FEM space)
        # (w,wt)+(grad u, grad wt) = <n, wt>, wt in H1
        # (grad w, grad ut)        = 0, ut in H01
        order = 2
        self.V = VectorH1(self.mesh, order=order)
        self.Q = VectorH1(self.mesh, order=order, dirichlet="bnd|b|r|t|l")
        self.X = FESpace([self.V,self.Q])
        self.sol = GridFunction(self.X)
        self.gInterface = GridFunction(self.X)
        self.gFixBnd = GridFunction(self.X)
        self.g0 = GridFunction(self.X)
        self.vector_lhs, self.vector_rhs = self.WeakForm()

    def WeakForm(self):
        # (w,wt)+(grad u, grad wt) = <n, wt>, wt in H1
        # (grad w, grad ut)        = 0, ut in H01
        w, u = self.X.TrialFunction()
        wt, ut = self.X.TestFunction()
        lhs = BilinearForm(self.X, symmetric=True)
        lhs += (InnerProduct(w,wt) + InnerProduct(grad(u),grad(wt)))*dx
        lhs += InnerProduct(grad(w),grad(ut))*dx
        rhs = LinearForm(self.X)
        rhs += -InnerProduct(grad(self.g0.components[1]),grad(wt))*dx + InnerProduct(specialcf.normal(2),wt)*ds
        return lhs, rhs

    def Solving(self,X_interface,dt):
        ## 设置Dirichlet边界条件：通过给定Interface上的值来求解BiHarmonic Extension
        ## 此时边界条件g0是[X,Q]上的函数，因此修改[1]component，位置延拓
        self.gInterface.components[1].Interpolate(X_interface,definedon=self.mesh.Boundaries('bnd'))
        self.gFixBnd.components[1].Interpolate(CoefficientFunction((x,y)),definedon=self.mesh.Boundaries('b|r|t|l'))
        self.g0.vec.data = self.gInterface.vec + self.gFixBnd.vec
        
        self.vector_lhs.Assemble()
        self.vector_rhs.Assemble()

        self.sol.vec.data = self.vector_lhs.mat.Inverse(self.X.FreeDofs()) * self.vector_rhs.vec + self.g0.vec
        self.VectorTmp.Interpolate(self.sol.components[1])
        self.Disp.vec.data = self.VectorTmp.vec - self.Initial_Pos.vec 

class X_DisLapExt(Extension):
    '''
        使用uniform weight的离散Laplace进行延拓，输入：
        mesh：网格
        Nfix：固定的网格点的数量，边界点的数量（其实可以从mesh中分析出来） 
        使用方法：
        * 初始化
        * Solving(Old_Coord(GridFunction))计算出displacement
    '''
    def __init__(self, mesh, Nfix=None, Ename='X') -> None:
        super().__init__(mesh, Ename)
        self.DisLap = self.Gen_Adjacent_Matrix()
        self.Nall = self.DisLap.shape[0]
        ## 边界点的数量
        if Nfix is not None:
            self.Nfix = Nfix
        else:
            bnd_index = []
            for el in self.mesh.Elements(BND):
                for v in el.vertices:
                    if v.nr not in bnd_index:
                        bnd_index.append(v.nr)
            self.Nfix = len(bnd_index)
        self.Lhs_Mat = self.Normal_Equation_Setting()
    
    def Gen_Adjacent_Matrix(self):
        u,v = self.fesscalar.TnT()
        lhs = BilinearForm(self.fesscalar)
        lhs += u*v*dx
        ## 通过邻接关系计算uniform weight下的discrete Laplace
        A = lhs.Assemble().mat
        A_data = list(A.COO())
        ## A_data[0]中每行每有一个非零元代表有一个相邻点
        Adjacent_nums = np.array(list(Counter(A_data[0]).values()))-1
        ## uniform weight，对角线为-1，其余均相等，总和为0
        N0 = len(A_data[0])
        A_data[2] = np.ones(N0)
        for ii in range(N0):
            row_num = A_data[0][ii]
            col_num = A_data[1][ii]
            if row_num != col_num:
                A_data[2][ii] = 1/Adjacent_nums[row_num]
            else:
                A_data[2][ii] = -1
        Lap_u = NGMO.myCOO(A_data[0],A_data[1],A_data[2],*(A.shape),'scipy')
        return Lap_u.copy()
            
    def Normal_Equation_Setting(self):
        ## 通过固定的节点用最小二乘依次计算其他节点的xy分量
        B = self.DisLap.tocsc()
        BlockDirichlet = NGMO.myCOO(list(range(self.Nfix)),list(range(self.Nfix)),list(np.ones(self.Nfix)),self.Nfix,B.shape[1],'scipy')
        BD = BlockDirichlet.tocsc()
        ## 稀疏的单位矩阵补零使得和B同样的列数
        Lhs_Mat = scipy.sparse.bmat([[B], [BD]])
        return Lhs_Mat

    def Solving(self,X_interface,dt=None):
        ## 分别利用法方程来求解两个分量：
        #   [Lu  ]   [new] = [ 0  ]
        #   [Im 0]   [ V ] = [Vbnd]
        Coords_interface = Pos_Transformer(X_interface)
        zero_vec = np.zeros(self.Nall).reshape(self.Nall,1)
        vx = scipy.sparse.bmat([[zero_vec],[Coords_interface[:,0].reshape(self.Nfix,1)]])
        Coords_new_x = LA.spsolve(self.Lhs_Mat.transpose()@self.Lhs_Mat,self.Lhs_Mat.transpose()@vx)
        vy = scipy.sparse.bmat([[zero_vec],[Coords_interface[:,1].reshape(self.Nfix,1)]])
        Coords_new_y = LA.spsolve(self.Lhs_Mat.transpose()@self.Lhs_Mat,self.Lhs_Mat.transpose()@vy)
        self.Disp.components[0].vec.data = BaseVector(Coords_new_x - self.Initial_Pos.components[0].vec.FV().NumPy())
        self.Disp.components[1].vec.data = BaseVector(Coords_new_y - self.Initial_Pos.components[1].vec.FV().NumPy())

def CircleSegs(index_set):
    tmp = index_set.reshape(len(index_set),1)
    tmp1 = tmp.copy()
    tmp1[:-1] = tmp[1:]
    tmp1[-1] = tmp[0]
    v = np.hstack((tmp,tmp1))
    return v

class HarmonicPullingBack(Extension):
    '''
        原始区域A，参数域Disk：B，两者有同样拓扑结构的网格，并且B的网格是非常好的，B0(在InitialParamMesh中生成)
        依据B的边界点给出A到B的边界的映射
        利用调和映射将A映射到B，得到B上的一个网格(这个过程需要A上的Deformed的mesh)，依据这个网格将B上的原始网格拉回到A上
        Coords_Origin_B: B上的原始网格
        Coords_Origin_A: A上的原始网格
        Coords_Updata_B: B上由A调和映射得到的网格
        mesh: 根据边界的移动需要改变的网格，和CoordMeshB具有相同的自由度排序
        mesh_A0: 求解
    '''
    
    def __init__(self,mesh):
        super().__init__(mesh,Ename="X")
        # self.UniformLapMod = X_DisLapExt(Mesh(mesh.ngmesh.Copy()))
        self.Bnd_Mapping, self.index = self.Initialize()
        self.bnd_num = len(self.index)
        self.CoordMeshB = self.SolvingHarmonicOnMesh()
        ## 每次更新的时候，根据节点坐标重新构造mesh_A0，特别是Delaunay三角化
        self.mesh_Mediate_A = Mesh(self.mesh.ngmesh.Copy())
        ## 然后将mesh_A0根据Bnd_Mapping映射到单位圆盘得到mesh_B0的网格
        self.mesh_B0 = Mesh(self.mesh.ngmesh.Copy())

    def Initialize(self):
        fesV = VectorH1(self.mesh,definedon=self.mesh.Boundaries("bnd"), order = 1)
        fesV = Compress(fesV)
        Pos_interface = GridFunction(fesV)
        Pos_interface.Interpolate(CoefficientFunction((x,y)),definedon=self.mesh.Boundaries(".*"))
        ## 边界点的坐标，第i行为编号为i的点的坐标
        Bnd_Coords = Pos_Transformer(Pos_interface)
        ## 取值在 -np.pi 到 np.pi 之间
        Bnd_theta = np.arctan2(Bnd_Coords[:,1],Bnd_Coords[:,0])
        ## index的第i个元素为从小往大排列第i个theta的index，例如第一个数为5，说明最小的是第5个点
        index = sorted(range(len(Bnd_theta)), key=lambda k: Bnd_theta[k])
        Circle_theta = np.linspace(-np.pi, np.pi, len(index)+1)[:-1]
        ## 构造边界映射
        Bnd_Mapping = np.zeros((len(index),2))
        for ii in range(len(index)):
            tmptheta = Circle_theta[ii]
            Bnd_Mapping[index[ii]] = np.array([np.cos(tmptheta),np.sin(tmptheta)])
        return Bnd_Mapping, index

    def SolvingHarmonicOnMesh(self,mesh=None):
        '''
            求解凸区域映射到disk上的Harmonic映射，凸区域的bnd的label为bnd
        '''        
        if mesh is not None:
            mymesh = mesh
        else: 
            mymesh = self.mesh
            print("Initialization Procedure!!")
        fesV = VectorH1(mymesh,definedon=mymesh.Boundaries("bnd"), order = 1)
        fesV = Compress(fesV)
        ## 存储边界映射的值，用来设置调和映射的Dirichlet边界条件
        Pos_interface_Disk = GridFunction(fesV)
        Pos_interface_Disk.vec.data = BaseVector(self.Bnd_Mapping.flatten('F'))
        fesVector = VectorH1(mymesh,order=1,dirichlet="bnd")
        g0 = GridFunction(fesVector)
        Harmonic_map = GridFunction(fesVector)
        ## 设置Dirichlet边界条件：基于A上的网格
        g0.Interpolate(Pos_interface_Disk,definedon=mymesh.Boundaries("bnd"))
        
        U,V = fesVector.TnT()
        lhs = BilinearForm(fesVector)
        lhs += InnerProduct(grad(U),grad(V))*dx
        rhs = LinearForm(fesVector)
        rhs += -InnerProduct(grad(g0),grad(V))*dx
        lhs.Assemble()
        rhs.Assemble()
        Harmonic_map.vec.data = g0.vec + lhs.mat.Inverse(inverse="pardiso",freedofs=fesVector.FreeDofs())*rhs.vec
        return Pos_Transformer(Harmonic_map)

    def PreMeshMod(self,Coords_Origin_A):
        '''
            根据节点坐标Coords_Origin_A来进行带限制边的Delaunay三角剖分生成mesh_Mediate_A
            Coords_Origin_A的节点排列顺序是边界上的节点排在最前面，将利用这些节点给出boundary segments
        '''
        data_set = {"vertices": Coords_Origin_A, "segments": CircleSegs(np.array(self.index))}
        tri = tr.triangulate(data_set,'p')
        ngmesh = ngm.Mesh(dim=2)
        pnums = []
        for v in Coords_Origin_A:
            pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(*v,0))))
        idx_dom = ngmesh.AddRegion("mat", dim=2)
        try:
            for tri_single in tri["triangles"]:
                ngmesh.Add(ngm.Element2D(idx_dom,[pnums[ii] for ii in tri_single]))
        except:
            print("{} is out of the total index of number {}".format(tri_single,Coords_Origin_A.shape[0]))
            np.save('error_data.npy',data_set,allow_pickle=True)
            print('index error occurred, error data are saved in error_data.npy!!')
            print(1/0)
        ## 添加边界
        idx_bnd = ngmesh.AddRegion("bnd", dim=1)
        try:
            for bnd in tri["segments"]:
                ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_bnd))
        except:
            print("{} is out of the total index of number {}".format(bnd,Coords_Origin_A.shape[0]))
            np.save('error_data.npy',data_set,allow_pickle=True)
            print('index error occurred, error data are saved in error_data.npy!!')
        self.mesh_Mediate_A = Mesh(ngmesh)

        ## 随mesh_A0变化的topological结构
        self.EL_VER_DICT = {}
        for el in self.mesh_Mediate_A.Elements(VOL):
            self.EL_VER_DICT[el.nr] = [v.nr for v in el.vertices]
        self.mesh_B0 = Mesh(self.mesh_Mediate_A.ngmesh.Copy())
        
        ## 通过求解mesh_Mediate_A到disk上的调和映射得到disk上的节点坐标
        Position = self.SolvingHarmonicOnMesh(mesh=self.mesh_Mediate_A)
        for i,vec in enumerate(self.mesh_B0.ngmesh.Points()):
            vec[0] = Position[i,0]
            vec[1] = Position[i,1]

        ## 构造Coords_Update_A
        ### 在新调整的mesh_B0上将self.CoordMeshB的节点值拉回到区域A上
        Coords_Update_A = np.zeros(Coords_Origin_A.shape)
        res = self.mesh_B0(self.CoordMeshB[:,0],self.CoordMeshB[:,1])
        for jj, it in enumerate(res):
            ## 单元index
            el_nr = it[-1]
            ## 单元的顺序节点index
            ver_nrs = self.EL_VER_DICT[el_nr]
            ## barycentric weight
            weight = [it[0],it[1],1-it[0]-it[1]]
            ## 拉回原始曲面上
            for ii in range(3):
                Coords_Update_A[jj] += weight[ii]*Coords_Origin_A[ver_nrs[ii]]
        Coords_Update_A[:self.bnd_num,:] = Coords_Origin_A[:self.bnd_num,:]
        self.Disp.vec.data = BaseVector(Coords_Update_A.flatten('F')) - self.Initial_Pos.vec

def SaveErrorData(err_data, errstr):
    '''
        错误的数据可以是ndarray格式，也可以是dict，其中每个value是ndarray格式
    '''
    print(errstr)
    np.save('error_data.npy',err_data,allow_pickle=True)

class HarmonicPullingBackQIR(Extension):
    def __init__(self,mesh,radi=0.5):
        super().__init__(mesh,Ename="X")
        self.radi = radi
        self.Bnd_Mapping, self.index, self.segs_outbnd, self.segs_inbnd = self.Initialize()
        self.bnd_index = self.segs_inbnd[:,0]
        self.bnd_num = len(self.index)
        ## 初始的CoordMeshB
        self.CoordMeshBIni = self.SolvingHarmonicOnMesh()
        ## 根据Jacobian调整的CoordMeshMod
        self.CoordMeshB = self.CoordMeshBIni.copy()
        ## 每次更新的时候，根据节点坐标重新构造mesh_A0，特别是Delaunay三角化
        self.mesh_Mediate_A = Mesh(self.mesh.ngmesh.Copy())
        ## 然后将mesh_A0根据Bnd_Mapping映射到单位圆盘得到mesh_B0的网格
        self.mesh_B0 = Mesh(self.mesh.ngmesh.Copy())

        ## 存储位置域上网格的拓扑连接关系，用以后文复原Laplace域中的网格用以存储
        self.Tri_Info = {}
        for el in self.mesh.Elements(VOL):
            self.Tri_Info[el.nr] = [v.nr for v in el.vertices]
        
        ## Get the boundary tag - segments dictionary from mesh
        bnd_dict = {}
        for el in self.mesh.Elements(BND):
            if el.mat in bnd_dict.keys():
                bnd_dict[el.mat].append([v.nr for v in el.vertices])
            else:
                bnd_dict[el.mat] = [[v.nr for v in el.vertices]]
        self.CoordMeshB_dict = {'vertices': self.CoordMeshB, 'triangles': self.Tri_Info, 'bnd_dict': bnd_dict}
        self.mesh_BMod = Tr2Mesh2d(self.CoordMeshB_dict)

    def Initialize(self):
        '''
            初始化从bnd到半径为radi的disk的边界的边界映射：
            将边界上的点按照弧度排序---arctan2，取值在 -np.pi 到 np.pi 之间
            然后按照这个顺序映射到以 self.radi 为半径的圆周的等距点上
        '''
        segments_set = {'b':[],'r':[],'t':[],'l':[],'bnd':[]}
        for el in self.mesh.Elements(BND):
            segments_set[el.mat].append([v.nr for v in el.vertices])
        for key,value in segments_set.items():
            segments_set[key] = np.array(value)
        segments_btrl = np.vstack((segments_set['b'],segments_set['t'],segments_set['r'],segments_set['l']))
        segments_bnd = segments_set['bnd']

        fesV = VectorH1(self.mesh,definedon=self.mesh.Boundaries("bnd"), order = 1)
        fesV = Compress(fesV)
        Pos_interface = GridFunction(fesV)
        Pos_interface.Interpolate(CoefficientFunction((x,y)),definedon=self.mesh.Boundaries("bnd"))
        ## 边界点的坐标，第i行为编号为i的点的坐标
        Bnd_Coords = Pos_Transformer(Pos_interface)
        ## 取值在 -np.pi 到 np.pi 之间
        Bnd_theta = np.arctan2(Bnd_Coords[:,1],Bnd_Coords[:,0])
        ## index的第i个元素为从小往大排列第i个theta的index，例如第一个数为5，说明最小的是第5个点
        index = sorted(range(len(Bnd_theta)), key=lambda k: Bnd_theta[k])
        Circle_theta = np.linspace(-np.pi, np.pi, len(index)+1)[:-1]
        ## 构造边界映射
        Bnd_Mapping = np.zeros((len(index),2))
        for ii in range(len(index)):
            tmptheta = Circle_theta[ii]
            Bnd_Mapping[index[ii]] = self.radi*np.array([np.cos(tmptheta),np.sin(tmptheta)])
        return Bnd_Mapping, index, segments_btrl, segments_bnd

    def SolvingHarmonicOnMesh(self,mesh=None):
        '''
            求解复合区域到矩形包含disk区域的Harmonic映射，
            mesh是个复合区域，边界有两部分构成，bnd以及矩形边界btrl，因此手动组装的mesh也需要指明这两类边界
        '''        
        if mesh is not None:
            mymesh = Mesh(mesh.ngmesh.Copy())
        else: 
            mymesh = Mesh(self.mesh.ngmesh.Copy())
            print("Initialization Procedure!!")
            mymesh.SetDeformation(self.mesh.deformation)
            
        fesVector = VectorH1(mymesh,order=1,dirichlet="bnd|b|r|t|l")
        g0 = GridFunction(fesVector)
        gInterface = GridFunction(fesVector)
        gFixBnd = GridFunction(fesVector)

        # gbnd是将mymesh的bnd边界映射到disk边界的分片线性函数，问题是如何按照正确的index构造
        # 首先构造限制定义在bnd上的有限元空间fesV，将其compress之后只剩下了bnd上的自由度，赋值，然后插值
        fesV = VectorH1(mymesh,definedon=mymesh.Boundaries("bnd"),order=1)
        fesV = Compress(fesV)
        ## 存储边界映射的值，用来设置调和映射的Dirichlet边界条件
        Pos_interface_Disk = GridFunction(fesV)
        Pos_interface_Disk.vec.data = BaseVector(self.Bnd_Mapping.flatten('F'))

        ## 设置Dirichlet边界条件：基于A上的网格
        gInterface.Interpolate(Pos_interface_Disk,definedon=mymesh.Boundaries("bnd"))
        gFixBnd.Interpolate(CF((x,y)),definedon=mymesh.Boundaries('b|t|r|l'))
        g0.vec.data = gInterface.vec + gFixBnd.vec

        Harmonic_map = GridFunction(fesVector)  
        U,V = fesVector.TnT()
        lhs = BilinearForm(fesVector)
        lhs += InnerProduct(grad(U),grad(V))*dx 
        rhs = LinearForm(fesVector) 
        rhs += -InnerProduct(grad(g0),grad(V))*dx   
        lhs.Assemble() 
        rhs.Assemble() 
        Harmonic_map.vec.data = g0.vec + lhs.mat.Inverse(inverse="pardiso",freedofs=fesVector.FreeDofs())*rhs.vec 
        return Pos_Transformer(Harmonic_map) 
    
    def CoordB_Recover(self):
        CoordsB = self.CoordMeshB
        Tri_Info = self.Tri_Info

        ngmesh = ngm.Mesh(dim=2)
        pnums = []
        for v in CoordsB:
            pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(*v,0))))
        idx_dom = ngmesh.AddRegion("mat", dim=2)
        for tri_single in Tri_Info.values():
            ngmesh.Add(ngm.Element2D(idx_dom,[pnums[ii] for ii in tri_single]))
        MeshB_Saved = Mesh(ngmesh)
        return MeshB_Saved

    def PreMeshMod(self,Coords_Origin_A):
        '''
            根据节点坐标Coords_Origin_A来进行带限制边的Delaunay三角剖分生成mesh_Mediate_A
            Coords_Origin_A的节点排列顺序是边界上的节点排在最前面，将利用这些节点给出boundary segments，既包含内边界，也包含外边界
        '''
        data_set = {"vertices": Coords_Origin_A, "segments": np.vstack((self.segs_outbnd,self.segs_inbnd))}
        tri = tr.triangulate(data_set,'-p')
        ngmesh = ngm.Mesh(dim=2)
        pnums = []
        for v in Coords_Origin_A:
            pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(*v,0))))
        idx_dom = ngmesh.AddRegion("mat", dim=2)
        for tri_single in tri["triangles"]:
            ngmesh.Add(ngm.Element2D(idx_dom,[pnums[ii] for ii in tri_single]))
        ## 添加内边界
        idx_bnd = ngmesh.AddRegion("bnd", dim=1)
        for bnd in self.segs_inbnd:
            ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_bnd))
        ## 添加外边界，简单起见均用label: 'b'
        idx_out_bnd = ngmesh.AddRegion("b", dim=1)
        for bnd in self.segs_outbnd:
            ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_out_bnd))

        self.mesh_Mediate_A = Mesh(ngmesh)

        ## 随mesh_A0变化的topological结构
        self.EL_VER_DICT = {}
        for el in self.mesh_Mediate_A.Elements(VOL):
            self.EL_VER_DICT[el.nr] = [v.nr for v in el.vertices]
        self.mesh_B0 = Mesh(self.mesh_Mediate_A.ngmesh.Copy())
        
        ## 通过求解mesh_Mediate_A到disk上的调和映射得到disk上的节点坐标
        Position = self.SolvingHarmonicOnMesh(mesh=self.mesh_Mediate_A)
        for i,vec in enumerate(self.mesh_B0.ngmesh.Points()):
            vec[0] = Position[i,0]
            vec[1] = Position[i,1]

        ## 从Mediate_A到mesh_B0的Jacobian
        fesV_MA = VectorH1(self.mesh_Mediate_A,order=1)
        Pos_B0 = GridFunction(fesV_MA)
        Pos_B0.vec.data = BaseVector(Position.flatten('F'))
        
        fes_MA = L2(self.mesh_Mediate_A,order=0)
        tmp_JacA2B = GridFunction(fes_MA)
        tmp_JacA2B.Interpolate(InnerProduct(grad(Pos_B0),grad(Pos_B0)))

        ## 构造Coords_Update_A
        ### 在新调整的mesh_B0上将self.CoordMeshB的节点值拉回到区域A上
        Coords_Update_A = np.zeros(Coords_Origin_A.shape)
        res = self.mesh_B0(self.CoordMeshB[:,0],self.CoordMeshB[:,1])
        for jj, it in enumerate(res):
            ## 单元index
            el_nr = it[-1]
            ## 单元的顺序节点index
            ver_nrs = self.EL_VER_DICT[el_nr]
            ## barycentric weight
            weight = [it[0],it[1],1-it[0]-it[1]]
            ## 拉回原始曲面上
            for ii in range(3):
                Coords_Update_A[jj] += weight[ii]*Coords_Origin_A[ver_nrs[ii]]
        Coords_Update_A[:self.bnd_num,:] = Coords_Origin_A[:self.bnd_num,:]
        self.Disp.vec.data = BaseVector(Coords_Update_A.flatten('F')) - self.Initial_Pos.vec

class HarmonicPullingBackQIRParametric():
    '''
        使用参数化的思想，将物理域的网格均通过参数域的均匀网格拉回得到，并且让物理域的网格与参数域的网格保持相同的拓扑连接关系：即作为后者的displacement存在即可。
        和之前的一大差别在于Physical Domain上的网格的节点编号继承圆周上的节点编号
    '''
    def __init__(self, msize, radi=0.5, rin=0.2, hInterface=None, Metric='Frob', No_Initial_Test=True):
        # No_Initial_Test：是否进行初始网格生成的测试，False，则将Harmonic拉回的disp和带扩散系数的Mod_Disp分开存储
        # Jacnorm: 'Frob' or 'Det'
        self.Diffusive_Metric_Tag = Metric
        self.Param_Radi = radi
        self.ParamMesh, ctr_point = CircleInRect(msize=msize, radi=radi, hInterface=hInterface)
        self.Coords_Param = np.array([v.point for v in self.ParamMesh.vertices])
        ## 获得参数域网格边界节点的index，用以后文强制映射到Quatrefoil上
        self.bnd_index = []
        for el in self.ParamMesh.Elements(BND):
            if el.mat == 'bnd':
                for v in el.vertices:
                    if v.nr not in self.bnd_index:
                        self.bnd_index.append(v.nr)
        self.bnd_num = len(self.bnd_index)
        self.Bnd_Mapping, self.segs_outbnd, self.segs_inbnd = self.Initialize()

        # 基于ParamMesh的有限元空间，因为Physical Mesh与Param具有相同的连接，也是Physical的有限元空间
        self.fesL0 = L2(self.ParamMesh, order=0, dirichlet='bnd|b|r|t|l')
        self.fesV = VectorH1(self.ParamMesh, order=1, dirichlet='bnd|b|r|t|l')
        self.Disp = GridFunction(self.fesV)
        self.tmp_JacA2B = GridFunction(self.fesL0)
        self.Initial_Area = Integrate(1,self.ParamMesh,element_wise=True).NumPy()
        
        self.PhysicalMesh, self.ctr_point = QuatrefoilInRect(msize=msize, rin=rin)
        ParamCoordUpdate = self.SolvingHarmonicOnMesh(self.PhysicalMesh,Init_opt=True)
        self.PhyCoordPullBack = CreatingMappingMesh(self.PhysicalMesh, ParamCoordUpdate, self.Coords_Param)
        ## 初始时刻将边界点映射回Quatrefoil interface上
        self.PhyCoordPullBack[self.bnd_index] = Pullback(0, 0, self.PhyCoordPullBack[self.bnd_index], Gen_Geometry_in_Splines3(self.ctr_point))
        # 初次的Disp，来自Harmonic拉回
        self.Disp.vec.data = BaseVector((self.PhyCoordPullBack-self.Coords_Param).flatten('F'))
        # 再次的Disp，来自Mod_By_Jac
        if No_Initial_Test:
            self.Disp.vec.data = BaseVector((Pos_Transformer(self.Mod_Phys_Mesh_By_Jac(self.Diffusive_Metric_Tag))-self.Coords_Param).flatten('F'))
        else:
            self.Mod_Disp = GridFunction(self.fesV)
            self.Mod_Disp.vec.data = BaseVector((Pos_Transformer(self.Mod_Phys_Mesh_By_Jac(self.Diffusive_Metric_Tag))-self.Coords_Param).flatten('F'))
        ## 此时的ParamMesh是有SetDeformation的，但不是最终到物理域网格的Disp

    def Initialize(self):
        segments_set = {'b':[],'r':[],'t':[],'l':[],'bnd':[]}
        for el in self.ParamMesh.Elements(BND):
            segments_set[el.mat].append([v.nr for v in el.vertices])
        for key,value in segments_set.items():
            segments_set[key] = np.array(value)
        segments_btrl = np.vstack((segments_set['b'],segments_set['t'],segments_set['r'],segments_set['l']))
        segments_bnd = segments_set['bnd']
        Bnd_Mapping = self.Coords_Param[np.sort(self.bnd_index)]
        return Bnd_Mapping, segments_btrl, segments_bnd

    def Mod_Phys_Mesh_By_Jac(self, Metric_Tag='Frob'):
        '''
            在参数域网格以及self.Disp的基础上再利用Jac映射进行调整：
            * 生成参数域网格上映射到deform之后网格的映射Pos_Update，生成Jacobian
            * 在PhysicalMesh（参数域网格+Disp）上求解带扩散系数的二阶方程 -- Dirichlet边界条件为保持PhysicalMesh的边界不变
        '''
        Coord_Update = Pos_Transformer(self.Disp) + self.Coords_Param
        ## PhysicalMesh的节点坐标以及Interface上的值（Dirichlet边界条件）
        Pos_Update = GridFunction(self.fesV)
        Pos_Update.vec.data = BaseVector(Coord_Update.flatten('F'))
        ## 在UnDeformed的 ParamMesh上生成将ParamMesh坐标映射到PhysicalMesh坐标的映射从而求Jacobian 
        # 再将标量函数Jacobian从ParamMesh提升到PhysicalMesh上（即SetDeformation）
        self.ParamMesh.UnsetDeformation()
        if Metric_Tag == 'Frob':
            # 扩散系数采用 <Jac, Jac> 其中Jac是Position映射的gradient（从参数域到物理域，因此必须先UnsetDeformation 
            self.tmp_JacA2B.Interpolate(InnerProduct(grad(Pos_Update),grad(Pos_Update)))
        ## 通过Deformation从ParamMesh得到PhysicalMesh 
        self.ParamMesh.SetDeformation(self.Disp)
        if Metric_Tag == 'Det':
            # 扩散系数采用 Det Jac，即Position映射对应的单元的面积变化率，参数域上单元的面积已经计算过，这里计算物理域上单元面积 
            Physical_Area = Integrate(1,self.ParamMesh,element_wise=True).NumPy()
            self.tmp_JacA2B.vec.data = BaseVector((Physical_Area/self.Initial_Area)**2)
        if Metric_Tag == 'TraceDet':
            A = grad(Pos_Update).trans*grad(Pos_Update)
            self.tmp_JacA2B.Interpolate(Trace(A)/Det(A))

        # 可以不用每次都重新生成
        g0 = GridFunction(self.fesV)
        U,V = self.fesV.TnT()
        lhs = BilinearForm(self.fesV)
        rhs = LinearForm(self.fesV)
        # lhs += InnerProduct(grad(U),grad(V))*(self.tmp_JacA2B)*dx
        # rhs += -InnerProduct(grad(g0),grad(V))*(self.tmp_JacA2B)*dx
        ## 
        lhs += InnerProduct(Inv(Id(2)+grad(Pos_Update).trans*grad(Pos_Update))*grad(U),grad(V))*dx
        rhs += -InnerProduct(Inv(Id(2)+grad(Pos_Update).trans*grad(Pos_Update))*grad(g0),grad(V))*dx
        sol = GridFunction(self.fesV)

        g0.Interpolate(CF((x,y)),definedon=self.ParamMesh.Boundaries('bnd|b|t|r|l'))
        lhs.Assemble()
        rhs.Assemble()
        sol.vec.data = lhs.mat.Inverse(self.fesV.FreeDofs(),inverse='pardiso')*rhs.vec + g0.vec
        return sol

    def SolvingHarmonicOnMesh(self,mesh,Init_opt=False):
        '''
            求解physical domain到 parametric domain的 Harmonic mapping
            mesh是个复合区域，边界有两部分构成，bnd以及矩形边界btrl，因此手动组装的mesh也需要指明这两类边界
            -----------
            interface条件的指定：
            * 初始时刻没有对应关系，因此将interface按照弧度顺序映射到圆周上，InterfaceCoord为nx2矩阵
            * 之后的physical domain作为 displacement，将interface映射到对应的圆周节点上
            输出Harmonic Mapping对应的节点坐标构成的nx2矩阵
        '''        
        # 第一种情形：输入mesh是初始的PhysicalMesh，求解Harmonic映到参数域的坐标，bnd_mat为 bnd|b|t|r|l
        # 第二种情形：输入mesh是新建的Mediate_A，求解Harmonic映射到参数域的坐标，bnd_mat为 bnd|b，两者都没有deformation
        PDmesh = Mesh(mesh.ngmesh.Copy())
        PDmesh.SetDeformation(mesh.deformation)
        
        PDfesV = VectorH1(PDmesh,order=1,dirichlet="bnd|b|r|t|l")
        g0 = GridFunction(PDfesV)
        gInterface = GridFunction(PDfesV)
        gFixBnd = GridFunction(PDfesV)

        # gbnd是将mymesh的bnd边界映射到disk边界的分片线性函数，问题是如何按照正确的index构造
        # 首先构造限制定义在bnd上的有限元空间fesV，将其compress之后只剩下了bnd上的自由度，赋值，然后插值
        PDfesV_bnd = Compress(VectorH1(PDmesh,definedon=PDmesh.Boundaries("bnd"),order=1))
        if Init_opt:
            # 生成初始网格到参数域的Laplace映射 --- 借助物理域中网格，因为此时物理域与参数域网格并不对应，手动设定interface的Dirichlet值
            # Physical domain的interface的坐标
            Pos_interface_Physic = GridFunction(PDfesV_bnd)
            Pos_interface_Physic.Interpolate(CF((x,y)), definedon=PDmesh.Boundaries('bnd'))
            self.Coord_Interface = Pos_Transformer(Pos_interface_Physic)
            # 按照theta顺序在圆周上等距排列
            theta_Interface = np.arctan2(self.Coord_Interface[:,1],self.Coord_Interface[:,0])
            index = sorted(range(len(theta_Interface)), key=lambda ii: theta_Interface[ii])
            theta_Circle = np.linspace(-np.pi,np.pi,len(theta_Interface)+1)[:-1]
            for ii in range(len(theta_Interface)):
                theta_Interface[index[ii]] = theta_Circle[ii]
            Bnd_Mapping = self.Param_Radi * np.array([[np.cos(theta), np.sin(theta)] for theta in theta_Interface])
        else:
            Bnd_Mapping = self.Bnd_Mapping
        ## 存储边界映射的值，用来设置调和映射的Dirichlet边界条件
        Pos_interface_Param = GridFunction(PDfesV_bnd)
        Pos_interface_Param.vec.data = BaseVector(Bnd_Mapping.flatten('F'))
        ## 设置Dirichlet边界条件：基于A上的网格
        gInterface.Interpolate(Pos_interface_Param,definedon=PDmesh.Boundaries("bnd"))
        gFixBnd.Interpolate(CF((x,y)),definedon=PDmesh.Boundaries('b|t|r|l'))
        g0.vec.data = gInterface.vec + gFixBnd.vec

        Harmonic_map = GridFunction(PDfesV)  
        U,V = PDfesV.TnT()
        lhs = BilinearForm(PDfesV)
        lhs += InnerProduct(grad(U),grad(V))*dx 
        rhs = LinearForm(PDfesV) 
        rhs += -InnerProduct(grad(g0),grad(V))*dx   
        lhs.Assemble() 
        rhs.Assemble() 
        Harmonic_map.vec.data = g0.vec + lhs.mat.Inverse(inverse="pardiso",freedofs=PDfesV.FreeDofs())*rhs.vec 
        return Pos_Transformer(Harmonic_map) 

    def PreMeshMod(self,Coords_Origin_A):
        '''
            根据节点坐标Coords_Origin_A来进行带限制边的Delaunay三角剖分生成mesh_Mediate_A
            Coords_Origin_A的节点排列顺序是边界上的节点排在最前面，将利用这些节点给出boundary segments，既包含内边界，也包含外边界
        '''
        data_set = {"vertices": Coords_Origin_A, "segments": np.vstack((self.segs_outbnd,self.segs_inbnd))}
        tri = tr.triangulate(data_set,'-p')
        ngmesh = ngm.Mesh(dim=2)
        pnums = []
        for v in Coords_Origin_A:
            pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(*v,0))))
        idx_dom = ngmesh.AddRegion("mat", dim=2)
        for tri_single in tri["triangles"]:
            ngmesh.Add(ngm.Element2D(idx_dom,[pnums[ii] for ii in tri_single]))
        ## 添加内边界
        idx_bnd = ngmesh.AddRegion("bnd", dim=1)
        for bnd in self.segs_inbnd:
            ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_bnd))
        ## 添加外边界，简单起见均用label: 'b'
        idx_out_bnd = ngmesh.AddRegion("b", dim=1)
        for bnd in self.segs_outbnd:
            ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_out_bnd))

        self.mesh_Mediate_A = Mesh(ngmesh)

        ## 随mesh_Mediate_A变化的topological结构
        self.EL_VER_DICT = {}
        for el in self.mesh_Mediate_A.Elements(VOL):
            self.EL_VER_DICT[el.nr] = [v.nr for v in el.vertices]
        self.mesh_B0 = Mesh(self.mesh_Mediate_A.ngmesh.Copy())
        
        ## 通过求解mesh_Mediate_A到disk上的调和映射得到disk上的节点坐标
        Position = self.SolvingHarmonicOnMesh(mesh=self.mesh_Mediate_A)
        for i,vec in enumerate(self.mesh_B0.ngmesh.Points()):
            vec[0] = Position[i,0]
            vec[1] = Position[i,1]

        ## 构造Coords_Update_A
        ### 在新调整的mesh_B0上将self.Coords_Param的节点值拉回到区域A上
        Coords_Update_A = np.zeros(Coords_Origin_A.shape)
        res = self.mesh_B0(self.Coords_Param[:,0],self.Coords_Param[:,1])
        for jj, it in enumerate(res):
            ## 单元index
            el_nr = it[-1]
            ## 单元的顺序节点index
            ver_nrs = self.EL_VER_DICT[el_nr]
            ## barycentric weight
            weight = [it[0],it[1],1-it[0]-it[1]]
            ## 拉回原始曲面上
            for ii in range(3):
                Coords_Update_A[jj] += weight[ii]*Coords_Origin_A[ver_nrs[ii]]
        Coords_Update_A[:self.bnd_num,:] = Coords_Origin_A[:self.bnd_num,:]
        self.Disp.vec.data = BaseVector((Coords_Update_A-self.Coords_Param).flatten('F'))
        # 再次拉回 -- 通过Mod_By_Jac
        self.Disp.vec.data = BaseVector((Pos_Transformer(self.Mod_Phys_Mesh_By_Jac())-self.Coords_Param).flatten('F'))

def CreatingMappingMesh(Physical_Mesh, Lap_Param_Pos, Coords_PullBack):
    '''
        Physical_Mesh 用来生成Laplace映前和映后的网格拓扑关系，并插值出网格映前节点位置
        Lap_Param_Pos 存储了Laplace映后的坐标，即网格映后节点位置
        Coords_PullBack 存储需要拉回的节点坐标
    '''
    # EL_VER_DICT 是Laplace映射域的拓扑关系，也是mesh_Rot自己的拓扑关系
    EL_VER_DICT = {}
    for el in Physical_Mesh.Elements(VOL):
        EL_VER_DICT[el.nr] = [v.nr for v in el.vertices]
    
    # mesh_Mapping 和求解完Laplace方程SetDeformation是一致的，不过作为mesh可以用来确定函数值
    mesh_Mapping = Mesh(Physical_Mesh.ngmesh.Copy())
    for i,vec in enumerate(mesh_Mapping.ngmesh.Points()):
        vec[0] = Lap_Param_Pos[i,0]
        vec[1] = Lap_Param_Pos[i,1]
        
    # mesh_Rot的初始位置
    fesV = VectorH1(Physical_Mesh,order=1)
    Pos_Ini = GridFunction(fesV)
    Pos_Ini.Interpolate(CF((x,y)))
    Coords_Origin = Pos_Transformer(Pos_Ini)
    
    # 利用mesh_Mapping将Coords_PullBack拉回并SetDeformation
    Coords_Update = np.zeros(Coords_PullBack.shape)
    res = mesh_Mapping(Coords_PullBack[:,0],Coords_PullBack[:,1])
    for jj, it in enumerate(res):
        ## 单元index
        el_nr = it[-1]
        ## 单元的顺序节点index
        ver_nrs = EL_VER_DICT[el_nr]
        ## barycentric weight
        weight = [it[0],it[1],1-it[0]-it[1]]
        ## 拉回原始曲面上
        for ii in range(3):
            Coords_Update[jj] += weight[ii]*Coords_Origin[ver_nrs[ii]]
#     Coords_Update[:self.bnd_num,:] = Coords_Origin[:self.bnd_num,:]
    return Coords_Update

def lift2sphere(mesh,R=2,Disp=None):
    '''
        输入ngsolve.comp.mesh（circle）将其lift到sphere上
    '''
    Coords = np.array([v.point for v in mesh.vertices])
    if not Disp is None:
        Coords += Pos_Transformer(Disp)
    ngmesh = ngm.Mesh(dim=3)
    pnums = []
    for v in Coords:
        z = np.sqrt(R**2-(np.linalg.norm(v))**2)
        if R == inf:
            # planar case
            z = 0
        pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(*v, z))))
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for el in mesh.Elements(VOL):
        el_v = [v.nr for v in el.vertices]
        ngmesh.Add(ngm.Element2D(idx_dom,[pnums[ii] for ii in el_v]))
    ### 设置Dirichlet边界条件
    idx_bnd = ngmesh.AddRegion("bnd", dim=1)
    for bnd in mesh.Elements(BND):
        bnd_v = [v.nr for v in bnd.vertices]
        ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd_v],index=idx_bnd))
    newmesh = Mesh(ngmesh)
    return newmesh

class BGNFixBnd():
    '''
        BGN method for mesh with nonempty boundary, but fixed.
        Thus, the finite element space is defined on a 2d surface in R3, with a fixed (Dirichlet) boundary of 1d
    '''
    def __init__(self, mesh):
        self.mesh = mesh
        self.fes = H1(self.mesh,order=1,definedon=self.mesh.Boundaries('.*'))
        self.fesV = VectorH1(self.mesh,order=1,definedon=self.mesh.Boundaries('.*'),dirichlet_bbnd='bnd')
        self.fesMix = self.fes*self.fesV
        self.Solution = GridFunction(self.fesMix)
        self.MethodName = 'BGN'
        ## 2d surface in 3d, mass lumping for triangles...
        ir = IntegrationRule(points = [(0,0), (1,0), (0,1)], weights = [1/6, 1/6, 1/6])
        self.ds_lumping = ds(intrules = { TRIG : ir })
        self.Position = GridFunction(self.fesV)
        self.Position.Interpolate(CoefficientFunction((x,y,z)),definedon=self.mesh.Boundaries('.*'))
        self.Disp = GridFunction(self.fesV)
        self.dim = 3
        self.vtk_bool = False
    
    def Initial_VTK(self,vtkname):
        T, Time_steps = 1,50   # useless params
        self.myvtk = Vtk_out_BND(T,min(50,Time_steps),vtkname)

    def WeakMCF(self):
        '''
            <D,chi n>^h = 0
                                        kappa n = -Lap X (球面的曲率是正的)
            <kappa n, eta>^h - <grad D, grad eta> = <grad Xold, grad eta>
            用来描述边界位移的D为Displacement，因为边界是fixed，所以D在边界上为0，恰好为homogeneous边界条件
        '''        
        kappa, D = self.fesMix.TrialFunction()
        chi, eta = self.fesMix.TestFunction()

        self.lhs = BilinearForm(self.fesMix)
        self.lhs += InnerProduct(D,specialcf.normal(3))*chi*self.ds_lumping
        
        self.lhs += InnerProduct(eta,specialcf.normal(3))*kappa*self.ds_lumping \
                    - InnerProduct(grad(D).Trace(),grad(eta).Trace())*ds

        self.rhs = LinearForm(self.fesMix)
        self.rhs += InnerProduct(grad(self.Position).Trace(),grad(eta).Trace())*ds

    def Solving(self,n):
        for ii in range(n):
            self.lhs.Assemble()
            self.rhs.Assemble()
            self.Solution.vec.data = self.lhs.mat.Inverse(self.fesMix.FreeDofs(),inverse="pardiso")*self.rhs.vec

            ## Solution.components[1] 代表的是X^m+1-X^m
            for i in range(self.dim):
                self.Disp.components[i].vec.data += self.Solution.components[1].components[i].vec.data  
                self.Position.components[i].vec.data += self.Solution.components[1].components[i].vec.data

            if self.vtk_bool:
                self.myvtk.Output(self.mesh,[],0,command='do')

            self.mesh.SetDeformation(self.Disp)

def Circle_Random_Mesh(N_Bnd=10,N_Int=40):
    '''
        单位圆中随机点通过Delaunay triangulation得到的三角剖分
    '''
    VER_BND = np.array([[np.cos(theta), np.sin(theta)] for theta in np.pi*np.arange(2*N_Bnd)/N_Bnd])
    r_Int = np.random.rand(N_Int)*0.99
    theta_Int = np.random.rand(N_Int)*2*np.pi
    VER_INT = np.array([[r*np.cos(theta), r*np.sin(theta)] for r, theta in zip(r_Int, theta_Int)])
    VER_ALL = np.vstack((VER_BND, VER_INT))
    data_set = {"vertices": VER_ALL, "segments": CircleSegs(np.arange(N_Bnd*2))}
    tri = tr.triangulate(data_set,'p')
    return tri

class Mod_Lift_Sphere():
    def __init__(self, mesh, R, disp):
        self.mesh = mesh
        self.fes = VectorH1(self.mesh, order=1)
        self.DeployDisp(R,disp)

    def DeployDisp(self,R,disp):
        '''
            if mesh has deformed, Object BGN_Mod should be reset.
        '''
        self.lift_mesh = lift2sphere(self.mesh, R=R, Disp=disp)
        self.BGN_Mod = BGNFixBnd(self.lift_mesh)
        self.BGN_Mod.WeakMCF()
        
    def Solving(self,Ntimes):
        self.BGN_Mod.Solving(Ntimes)
    
    def ProjDisp(self):
        Proj_Disp = Pos_Transformer(self.BGN_Mod.Disp,3)
        Proj_Disp = Proj_Disp[:,:2]
        return Proj_Disp

def Tr2Mesh2d(Tr_Info_Dict):
    VER_ALL = Tr_Info_Dict['vertices']
    ngmesh = ngm.Mesh(dim=2)
    pnums = []
    for v in VER_ALL:
        pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(*v,0))))
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for tri_single in Tr_Info_Dict["triangles"].values():
        ngmesh.Add(ngm.Element2D(idx_dom,[pnums[ii] for ii in tri_single]))
    ## 添加边界
    if 'bnd_dict' in Tr_Info_Dict.keys():
        for key, value in Tr_Info_Dict['bnd_dict'].items():
            idx_bnd = ngmesh.AddRegion(key, dim=1)    
            for bnd in value:
                ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_bnd))
    elif 'segments' in Tr_Info_Dict.keys():
        # legacy
        idx_bnd = ngmesh.AddRegion("bnd", dim=1)
        for bnd in Tr_Info_Dict['segments']:
            ngmesh.Add(ngm.Element1D([pnums[ii] for ii in bnd],index=idx_bnd))
    mesh_2d = Mesh(ngmesh)
    return mesh_2d