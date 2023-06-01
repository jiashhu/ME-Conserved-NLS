from netgen.meshing import Mesh, MeshPoint, Element1D, FaceDescriptor, Element0D, Element2D
from netgen.csg import Pnt
import numpy as np
from netgen.geom2d import *

class Mesh1d():
    def __init__(self,lp,rp,N,periodic=False) -> None:
        assert(lp<rp)
        self.points_set = np.linspace(lp,rp,N+1)
        self.h = (rp-lp)/N
        self.ngmesh = Mesh(dim=1)
        self.N = N
        self.periodic = periodic
        self.GenerateMesh()

    def GenerateMesh(self):
        pnums = []
        for i in range(self.N + 1):
            pnums.append(self.ngmesh.Add(MeshPoint(Pnt(self.points_set[i],0,0))))
        self.ngmesh.SetMaterial(1, "mat")
        for j in range(self.N):
            # or using Element1D([pnums[j],pnums[j+1]],index=1)
            self.ngmesh.Add(Element1D([pnums[j], pnums[j + 1]],index=1))
        self.ngmesh.Add (Element0D(pnums[0], index=1))
        self.ngmesh.Add (Element0D(pnums[-1], index=2))
        if self.periodic:
            self.ngmesh.AddPointIdentification(pnums[0],pnums[-1],1,2)

class PeriodicSquare():
    '''
        Omega = [0,1]x[0,1]
    '''
    def __init__(self,h,L=1/2,periodic=False):
        self.h = h
        self.L = L
        self.periodic = periodic
        self.GenerateMesh()

    def GenerateMesh(self):
        periodic = SplineGeometry()
        pnts = [ (-self.L,-self.L), (self.L,-self.L), (self.L,self.L), (-self.L,self.L) ]
        pnums = [periodic.AppendPoint(*p) for p in pnts]

        ldown = periodic.Append ( ["line", pnums[0], pnums[1]],bc="outer")
        # This should be our master edge so we need to save its number.
        lright = periodic.Append ( ["line", pnums[1], pnums[2]], bc="periodic")
        periodic.Append ( ["line", pnums[3], pnums[2]], leftdomain=0, rightdomain=1, copy=ldown, bc="outer")
        # Minion boundaries must be defined in the same direction as master ones,
        # this is why the the point numbers of this spline are defined in the reverse direction,
        # leftdomain and rightdomain must therefore be switched as well!
        # We use the master number as the copy argument to create a slave edge.
        periodic.Append ( ["line", pnums[0], pnums[3]], leftdomain=0, rightdomain=1, copy=lright, bc="periodic")
        self.ngmesh = periodic.GenerateMesh(maxh=self.h)

class PeriodicSquareUniform():
    '''
        Omega = [-L,L]x[-L,L]; Rot_opt: rotation pi/4 clockwise
    '''
    def __init__(self, N, L=1/2, periodic=False, Rot_opt=False):
        self.N = N
        self.L = L
        self.periodic = periodic
        self.Rot_opt = Rot_opt
        self.GenerateMesh()
    
    def GenerateMesh(self):
        self.ngmesh = Mesh()
        self.ngmesh.dim = 2

        pnums = []
        for i in range(self.N + 1):
            for j in range(self.N + 1):
                xij = -self.L + i*2*self.L/self.N
                yij = -self.L + j*2*self.L/self.N
                if self.Rot_opt:
                    rot_x = np.sqrt(2)/2*(xij+yij)
                    rot_y = np.sqrt(2)/2*(-xij+yij)
                    xij, yij = rot_x, rot_y
                pnums.append(self.ngmesh.Add(MeshPoint(Pnt(xij, yij, 0))))
                if self.periodic:
                    if j == 0:
                        jslave = pnums[-1]
                    if j == self.N:
                        self.ngmesh.AddPointIdentification(pnums[-1],jslave,identnr=1,type=2)
        
        if self.periodic:
            for i in range(self.N+1):
                islave = pnums[i]
                self.ngmesh.AddPointIdentification(pnums[-(self.N+1)+i],islave,identnr=2,type=2)

        idx = self.ngmesh.AddRegion("mat",dim=2)
        for j in range(self.N):
            for i in range(self.N):
                self.ngmesh.Add(Element2D(idx, [pnums[i + j * (self.N + 1)], pnums[i + (j + 1) * (self.N + 1)], pnums[i + 1 + j * (self.N + 1)]]))
                self.ngmesh.Add(Element2D(idx, [pnums[i + (j + 1) * (self.N + 1)], pnums[i + 1 + (j + 1) * (self.N + 1)], pnums[i + 1 + j * (self.N + 1)]]))