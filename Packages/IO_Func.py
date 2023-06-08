import os
import numpy as np
from ngsolve import *
import ngsolve
from pyevtk.hl import polyLinesToVTK
from netgen.meshing import Mesh, MeshPoint, Element1D, FaceDescriptor, Element0D, Element2D
from netgen.csg import Pnt
import pytz
import datetime
try:
    import matplotlib as mpl
except:
    print('no matplotlib is found')
from cycler import cycler

def Get_File_List(file_path):
    '''
        output: list of file name, ranked by created time
    '''
    dir_list = os.listdir(file_path)
    if not dir_list:
        print('File path is empty!!')
        return None
    else:
        dir_list = sorted(dir_list, key=lambda x: os.path.getctime(os.path.join(file_path,x)))
        return dir_list

def Transformh(L_Total,Nset):
    h = [L_Total/Ni for Ni in Nset]
    return h
    
def GenerateLim(minh,maxh,minE,maxE,adj_param,opt='centeral'):
    if opt == 'max':
        ylim1,ylim2 = 10**np.floor(np.log(minE)/np.log(10)),10**np.ceil(np.log(maxE)/np.log(10))
        ratio = ylim2/ylim1
        xlim1,xlim2 = np.sqrt(minh*maxh)/np.sqrt(ratio), np.sqrt(minh*maxh)*np.sqrt(ratio)
    elif opt == 'centeral':
        cx = np.log(np.sqrt(minh*maxh))/np.log(10)
        cy = np.log(np.sqrt(minE*maxE))/np.log(10)
        hx = np.log(maxh/minh)/np.log(10)
        hy = np.log(maxE/minE)/np.log(10)
        ylim1, ylim2 = 10**(cy-hy/2), 10**(cy+hy/2)
        xlim1, xlim2 = 10**(cx-hx/2), 10**(cx+hx/2)
    return xlim1,xlim2,ylim1,ylim2

def AnalyzeParam(h,err,order,adj_param,L_Total,xyCoef=None,add_order=0,opt='centeral'):
    # big h corresponds to big error, adj_param<1
    minE,maxE = min(err),max(err)
    minh,maxh = min(h), max(h)
    if np.log(err[0]/err[-1])/np.log(h[0]/h[-1]) > 0:
        pass
    else:
        h = Transformh(L_Total,h)
        minh,maxh = min(h), max(h)
    if xyCoef is not None:
        xlim1,xlim2,ylim1,ylim2 = GenerateLim(minh,maxh,minE,maxE,adj_param,opt=opt)
        xlim1,xlim2,ylim1,ylim2 = np.array([xlim1,xlim2,ylim1,ylim2])*xyCoef
    else:
        xlim1,xlim2,ylim1,ylim2 = GenerateLim(minh,maxh,minE,maxE,adj_param,opt=opt)

    hstr = [format(hi,'0.2e') for hi in h]
    refh = h.copy()
    refE = [maxE*adj_param**(float(order)+add_order)*(hi/maxh)**(float(order)+add_order) for hi in h]
    return h, err, refh, refE, hstr, [xlim1,xlim2], [ylim1,ylim2]

class CanonicalErrPlot():
    def __init__(self,Res,Zoom=None,XYTick=None,Legend=None,adj_p=0.8,opt='centeral') -> None:
        self.Res = Res
        self.index = 1
        self.adj_param = adj_p
        self.Zoom = Zoom
        self.XYTick = XYTick
        self.Legend = Legend
        self.add_order = self.Res['Params']['add_order']
        self.opt = opt
        self.var = self.Res['Params']['var']

    def ReadParam(self):
        self.Err_Data = self.Res['Converg']
        self.fig_name, self.Param, self.Err = [self.Res['Params'][index] for index in ['Name','Param','Err']]
        self.n_order = len(self.Err_Data)
        if type(self.adj_param) is float:
            self.adj_param = [self.adj_param]*self.n_order
        self.L_Tot = self.Res['Params']['Param_Total']

    def DrawParam(self,fig):
        set_ylabel = True
        for key,val in self.Err_Data.items():
            order = key.split('_Conv_o')[-1]
            ax = fig.add_subplot(1,self.n_order,self.index)
            h, err = val[self.Param], val[self.Err]
            # self.index starts from 1
            h,err,refh,refE,hstr,xlim,ylim = AnalyzeParam(h,err,order,self.adj_param[self.index-1],self.L_Tot,add_order=self.add_order,opt=self.opt)
            ax.loglog(h,err,label=self.Legend[order][0])
            ax.loglog(refh,refE,label=self.Legend[order][1])
            ax.grid(True,which="both",ls="--")            
            try:
                if self.Zoom is not None:
                    xlim = np.array(xlim)*np.array(self.Zoom[order]['x'])
                    ylim = np.array(ylim)*np.array(self.Zoom[order]['y'])
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
            except:
                pass
            ax.set_xlabel(self.Res['Params']['xlabel'])
            if set_ylabel:
                ax.set_ylabel('$H^1 error$')
                set_ylabel = False
            ax.set_title('='.join([self.var,order]))

            try:
                if self.XYTick[order] is not None:
                    ax.set_xticks(self.XYTick[order]['xtick'])
                    ax.set_xticklabels(self.XYTick[order]['xticklabels'])
                    ax.set_yticks(self.XYTick[order]['ytick'])
                    ax.set_yticklabels(self.XYTick[order]['yticklabels'])
                    ax.set_yticklabels([],minor=True)
            except:
                ax.set_xticklabels([],minor=True)
                ax.set_yticklabels([],minor=True)
            ax.legend(loc='lower right')
            self.index += 1

class PlotLineStyleObj():
    def __init__(self) -> None:
        self.ColorCycle = None
        self.MarkerCycle = None
        self.LinestyleCycle = None

    def SetColorCycle(self):
        self.ColorCycle = cycler('color',[
                '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', 
                '#bcbd22', '#17becf', 'k', 'b'
            ])

    def SetMarkerCycle(self):
        self.MarkerCycle = cycler('marker',[
                '.',',','o','v'
            ]*3)
    
    def SetLinestyleCycle(self):
        self.LinestyleCycle = cycler('linestyle',['-','--','-.',':']*3)

    def SetPropCycle(self, linewidth=2, markersize=16,
        markeredgewidth = 2, fontsize=14):
        mpl.rcParams['axes.prop_cycle'] = self.ColorCycle
        for mycycler in [self.LinestyleCycle,self.MarkerCycle]:
            if mycycler is not None:
                mpl.rcParams['axes.prop_cycle'] +=  mycycler
        mpl.rcParams['lines.markersize'] = markersize
        mpl.rcParams['lines.linewidth'] = linewidth
        mpl.rcParams['lines.markeredgewidth'] = markeredgewidth
        mpl.rcParams['lines.markerfacecolor'] = 'none'
        mpl.rcParams.update({'font.size': fontsize})

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

def PVD_Generate(pvd_path,folder_path_set:list,pvd_name,T_end_set:list=[np.inf]):
    '''
        file strcture:
        - pvd_path (single path for pvd file)
            - folder1 (for series of vtk files)
            - folder2
            ...

        the subfolders contain data from different period, and folder_path_set collects all these paths.

        If there are several subfolders, choose vtus in each one according to T_end, and collect in the pvd file. 
        In the folder, time infomation is included in Rel_Mapping.npy (dictionary: vtu_name: time)
    '''
    assert len(T_end_set) == len(folder_path_set)
    T_begin = 0
    pvd_file_list_begin = [
            '<?xml version="1.0"?>',
            '<VTKFile type="Collection" version="0.1"',
            '    byte_order="LittleEndian"',
            '    compressor="vtkZLibDataCompressor">',
            '<Collection>']
    pvd_file_list_content = []
    for vtu_folder_rel_path, T_end in zip(folder_path_set,T_end_set):
        abs_path = os.path.join(pvd_path,vtu_folder_rel_path)
        # to match time steps with vtk number, save this pair in pvd file
        try:
            Map_Data = np.load(os.path.join(abs_path,'Rel_Mapping.npy'),allow_pickle=True).item()
            time_set = np.array(list(Map_Data.values()))
            vtu_name = list(Map_Data.keys())
        except:
            # Arrange vtu files in alphabetical order by file name
            print('No Mapping data!')
            vtu_name = [fname.split('.')[0] for fname in Get_File_List(abs_path) if fname.endswith('vtu')]
            time_set = np.array(list(range(len(vtu_name))))
        # find the index of time steps in the range (default time_set increases)
        index = (time_set>=T_begin) & (time_set<T_end)
        time_period = time_set[index]
        vtu_period = [v_i for v_i, b_v in zip(vtu_name,index) if b_v == True]
        T_begin = T_end

        for vtk_file, t_vtk in zip(vtu_period,time_period):
            mystr1 = '    <DataSet timestep="{}" group="" part="0"'.format(t_vtk)
            mystr2 = '            file="{}.vtu"/>'.format(os.path.join(vtu_folder_rel_path,vtk_file))
            pvd_file_list_content.append(mystr1)
            pvd_file_list_content.append(mystr2)
            
    pvd_file_list_end = [
            '</Collection>',
            '</VTKFile>'
        ]    
    pvd_file_all = pvd_file_list_begin + pvd_file_list_content + pvd_file_list_end
    with open(os.path.join(pvd_path,pvd_name), 'w') as f:
        for single_line in pvd_file_all:
            f.write(single_line + '\n')

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
    
    def GenerateMapping(self, filepath, LoadT=0):
        '''
            LoadT = 0 version
        '''
        vtu_list = [fname for fname in Get_File_List(filepath) if fname.endswith('vtu')]
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
        PVD_Generate(pvd_path=self.pathname,folder_path_set=[''], pvd_name=pvd_name)

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
            Input: vertices of 1d curve (in anticlockwise order)
            Output: vtk file
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
        self.mesh     = ngsolve.Mesh(Mesh_Obj.ngmesh)
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
        
def LogTime(timezone='Asia/Shanghai'):
    format = "%Y-%m-%d %H:%M:%S %z"
    a_datetime = datetime.datetime.now(pytz.timezone(timezone))
    datetime_string = a_datetime.strftime(format) + " " + timezone
    return datetime_string

class Printer():
    '''
        myPrinter = Printer(10)
        myPrinter.PrintDuring(told,T,T_begin=0) 
    '''
    def __init__(self,n) -> None:
        self.nprint = n
        self.n_now = 0

    def PrintDuring(self,told,T,T_begin=0):
        if told>=T_begin+(T-T_begin)*self.n_now/self.nprint:
            print("{}-Finished {} per cent".format(LogTime(),self.n_now/self.nprint*100))
            self.n_now += 1

def PosRemove0(mystr:str):
    delet_terms = ['e+00']
    for terms in delet_terms:
        if terms in mystr:
            mystr = mystr.split(terms)[0]
    stripzero = True
    if 'e+' in mystr:
        keyw = 'e+'
    elif 'e-' in mystr:
        keyw = 'e-'
    else:
        stripzero = False
    if stripzero:
        mystr = keyw.join(res.strip('0') for res in mystr.split(keyw))
    return mystr

def GenerateTick_Label(tick_val:list,opt):
    tick_name = opt+'tick'
    tlabel_name = tick_name+'labels'
    tlabel = [PosRemove0(format(v,'0.0e')) for v in tick_val]
    mydict = dict(zip([tick_name,tlabel_name],[tick_val,tlabel]))
    return mydict

def GenerateXYTick(val:list):
    Res = []
    for x,y in val:
        res = GenerateTick_Label(x,'x')
        res.update(GenerateTick_Label(y,'y'))
        Res.append(res)
    return Res