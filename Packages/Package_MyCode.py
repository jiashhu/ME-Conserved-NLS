import os
import numpy as np
import math
import re
try:
    import pytz
    from cycler import cycler
    import logging
    from PIL import Image
    import matplotlib as mpl
    from line_profiler import LineProfiler
except:
    pass
import datetime

class SliceableDict(dict):
    '''
        Initialize the same as the normal dictionary: SliceableDict(Param_dict) 
        by sequence = {key:val for key,val in zip(range(5),range(2,6))}
    '''
    def slice(self, *keys):
        return {k: self[k] for k in keys}

def PNG2PDF(BaseDir):
    flist = FO.Get_File_List(BaseDir)
    for fname in flist:
        if fname.endswith('.png'):
            newname = '-'.join(fname.split('.')[:-1])
            rgba = Image.open(os.path.join(BaseDir,fname))
            rgb = Image.new('RGB', rgba.size, (255, 255, 255))  # white background
            rgb.paste(rgba, mask=rgba.split()[3])               # paste using alpha channel as mask
            rgb.save(os.path.join(BaseDir,'{}.pdf'.format(newname)), "PDF",
                                resolution=100.0, save_all=True)

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
            assert key.startswith('order')
            order = key.split('order')[-1]
            assert all([kws in val.keys() for kws in [self.Param,self.Err]])
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
    
def myfloat(s):
    try:
        tmp = float(s)
        return tmp
    except:
        return s

def LogSet(fname = 'example.log'):
    logger = logging.getLogger()
    logger.setLevel('DEBUG')
    BASIC_FORMAT = "%(asctime)s:%(levelname)s:%(message)s"
    DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(BASIC_FORMAT, DATE_FORMAT)
    chlr = logging.StreamHandler() # 输出到控制台的handler
    chlr.setFormatter(formatter)
    chlr.setLevel('INFO')  # 也可以不设置，不设置就默认用logger的level
    fhlr = logging.FileHandler(fname) # 输出到文件的handler
    fhlr.setFormatter(formatter)
    logger.addHandler(chlr)
    logger.addHandler(fhlr)
    return logger

def LogVariables(dict,Prompt='Params',NumPre=False,**kwargs):
    for key_add in kwargs:
        dict[key_add] = kwargs[key_add]
    mykey = list(dict.keys())
    if NumPre:
        output = ','.join([key+'='+format(val,'0.6e') for key,val in dict.items()])
    else:
        output = ','.join([key+'='+str(dict[key]) for key in mykey])
    return '{}: {}'.format(Prompt,output)

class LogParse():
    '''
        Parse the log info containing the keyword: default is = 'Initialize'
    '''
    def __init__(self, fname = 'example.log', keyword='Initialize') -> None:
        self.fname = fname
        self.Parameter_done = []
        with open(fname,'rb') as f:
            myline = 'begin'
            while myline:
                myline = f.readline().decode().replace('\n','').replace(' ','')
                ## Parameter format: keyword .... : key1=val1, ... keyn=valn
                if keyword in myline:
                    data = re.split(':',myline)[-1]
                    res = re.split('=|,',data)
                    res = [myfloat(num) for num in res]
                    param_done = {key:val for key,val in zip(res[::2],res[1::2])}
                    self.Parameter_done.append(param_done)
    
    def check(self,checkdict):
        check_res = False
        for Params in self.Parameter_done:
            if checkdict.items() <= Params.items():
                print('{} already done.'.format(', '.join(['='.join([key,str(val)]) 
                                        for key,val in checkdict.items()])))
                check_res = True
                break
        return check_res

def do_profile(follow=None):
    """
    使用line_profiler创建性能分析装饰器
    follow列表选择要追踪的函数，如果为空，则全部分析
    用例：
    def num_range(n):
        for x in range(n):
            yield x

    @do_profile(follow=[num_range])
    def expensive_function():
        for x in num_range(1000):
            _ = x ** x
        return 'OK!'

    result = expensive_function()
    """
    if follow is None:
        follow = list()
    def inner(func):
        def profiled_func(*args, **kwargs):
            profiler = LineProfiler()
            try:
                profiler.add_function(func)
                for f in follow:
                    profiler.add_function(f)
                profiler.enable_by_count()
                return func(*args, **kwargs)
            finally:
                profiler.print_stats()
        return profiled_func
    return inner

def LogTime(timezone='Asia/Shanghai'):
    format = "%Y-%m-%d %H:%M:%S %z"
    a_datetime = datetime.datetime.now(pytz.timezone(timezone))
    datetime_string = a_datetime.strftime(format) + " " + timezone
    return datetime_string

class FO():
    def Edit_Name(FoldPath,FileList,Prefix='',DelItem=''):
        '''删除文件名中DelItem, 增加前缀Prefix'''
        for fname in FileList:
            NewName = fname
            ModFlag = 0
            if DelItem:
                if DelItem in NewName:
                    NewName = NewName.replace(DelItem,'')
                    ModFlag = 1
            if Prefix:
                if not NewName.startswith(Prefix):
                    NewName = ''.join([Prefix,NewName])
                    ModFlag = 1
            if ModFlag:
                os.rename(os.path.join(FoldPath,fname),os.path.join(FoldPath,NewName))

    def Create_Data_Fold(DF_Path_Rel,MeshDirName,MName):
        '''
            建立数据存储文件夹: DF_Path_Rel + Mesh名 + 方法名，
            ./Data_MCF_Split/Dumbbell_Singular/Dziuk，
            并将初始网格(vol格式) copy进去
        '''
        Datafilepath = '{}/{}/{}/'.format(DF_Path_Rel,MeshDirName,MName)
        if not os.path.exists(Datafilepath):
            os.makedirs(Datafilepath)
            print('目标文件夹：{}不存在，但是已经创建'.format(Datafilepath))
        else:
            print('目标文件夹：{}已存在'.format(Datafilepath))
        # 将当前文件夹下的Mesh文件复制到Data文件夹下
        MeshName = MeshDirName+'.vol'
        CopyMName = Datafilepath+MeshName
        if not os.path.exists(CopyMName):
            try:
                assert(os.path.exists('./'+MeshName))
                os.system('cp {} {}'.format(MeshName,CopyMName))
                print('将初始网格copy')
            except:
                print('当前文件夹下没有初始mesh文件')
        else: 
            print('初始网格已经在文件夹中')
        print('DataFold is {}'.format(Datafilepath))
        return Datafilepath

    def Get_File_List(file_path):
        '''
            获取当前文件夹下file
            输入: 文件路径
            输出: list of file name (按创建时间排序)   
        '''
        dir_list = os.listdir(file_path)
        if not dir_list:
            print('File path is empty!!')
            return None
        else:
            dir_list = sorted(dir_list, key=lambda x: os.path.getctime(os.path.join(file_path,x)))
            return dir_list

    def DictValueCheck(Dict):
        '''
            将可执行的value转化为执行后的结果
        '''
        for key, val in Dict.items():
            try:
                Dict[key] = eval(val)
            except:
                pass
        return Dict

    def Generate_Name(ParamDict,MName,n=4):
        '''
        输入参数: 参数Dict, 参数名少于2个字母, 并将MName词条排到第一个--值需要为str,
        参数值如果是数字, 采用科学计数表示, n表示小数点后的位数
        '''
        Param_list = [ParamDict[MName]]
        for key, val in ParamDict.items():
            if key == MName:
                continue
            else:
                if type(val) == str:
                    pass
                ## 处理万一不是str类型的值的情况
                elif type(val) == int:
                    val = str(val)
                elif type(val) != str:
                    val = MO.Num2Sci(val,n)
                Param_list.append(key+'-'+val)
        res = '-'.join(Param_list).replace('.','_').replace('/',':')
        print("File Name is {}".format(res))
        return res

    def Parse_Name(DictPath,label='T1',suffix='.gz'):
        '''解析DictPath下数据文件名中的label, 返回Dict'''
        MeshInfoDict = {}
        info = []
        for fname in FO.Get_File_List(DictPath):
            if fname.endswith(suffix):
                info = fname.split('.')[0].split('-')[1:]
                if label in info:
                    ind = info.index(label)
                    MeshInfoDict[fname] = info[ind+1]
        return MeshInfoDict

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
                vtu_name = [fname.split('.')[0] for fname in FO.Get_File_List(abs_path) if fname.endswith('vtu')]
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

    def Copy_Code(OriginDir, TargetDir, CodeName, suffix=''):
        target_path = os.path.join(TargetDir, CodeName)
        present_path = os.path.join(OriginDir, CodeName)
        copy_opt = False
        if not OriginDir==TargetDir:    
            if os.path.exists(target_path):
                val = input('Check file, override (yes or no):')
                if val == 'yes':
                    copy_opt = True
                    os.system('rm {}'.format(target_path))
                else:
                    print('The code file is not copied!')
            else:
                copy_opt = True
            if copy_opt:
                codename, codesuff = CodeName.split('.')
                os.system('cp {} {}'.format(present_path,os.path.join(TargetDir,codename+suffix+'.'+codesuff)))
                print('Backup successfully!')
        else:
            print('The code is in the present directory.')
class MO():
    def Num2Sci(x,n):
        '''
        将numeric x保存小数点后n位有效数字的科学计数法输出成文本
        '''
        return format(x,'0.{}e'.format(n)).replace('.','_').replace('-','_')

class Pic():
    def Get_Boundary(img):
        '''输入Image.Open的结果，输出有效部分的边界'''
        im = np.array(img.convert('L'))
        grad = np.gradient(im)
        gradx,grady = grad
        ## 按列求和（求行和），找出第一个和最后一个含有非零元的行
        res = np.where(abs(gradx).sum(axis=1)!=0)[0]
        # 左上角为原点，往下为x轴，往右为y轴
        x0,x1 = min(res),max(res)
        res = np.where(abs(grady).sum(axis=0)!=0)[0]
        y0,y1 = min(res),max(res)
        return x0,x1,y0,y1

    def Cut(BaseFold,RelPath,DirectFold,auto=True,w=(0.8,0.8),h=(0.6,0.6)):
        '''
        截取BaseFold+RelPath中的图片，保存到DirectFold中，
        默认是自动截取灰度非零的部分，
        width, height 为长与宽的比例，
        eps帮助见 https://stackoverflow.com/questions/47398291/saving-to-eps-not-supported-in-python-pillow
        '''
        if not DirectFold.endswith('/'):
            DirectFold += '/'
        for filename in RelPath:
            img = Image.open(BaseFold+filename)
            ## Eps设置
            if img.mode in ('RGBA', 'LA'):
                print('Current figure mode "{}" cannot be directly \
                        saved to .eps and should be converted \
                        (e.g. to "RGB")'.format(img.mode))
                img = img.convert('RGB')

            OutName = DirectFold+filename.split('.')[0]+'.eps'
            if os.path.exists(OutName):
                print("文件：{}已经存在".format(OutName))
                continue
            else:
                if auto:
                    upper, lower, left, right = Pic.Get_Boundary(img)
                else:
                    left = (1-w[0])/2*img.size[0]
                    right = (1+w[1])/2*img.size[0]
                    upper = (1-h[0])/2*img.size[1]
                    lower = (1+h[1])/2*img.size[1]
                cropped = img.crop((left, upper, right, lower))  # (left, upper, right, lower)
                cropped.save(OutName)
                print(OutName)
                img.close()
                cropped.close()
    
    def Resize_Canvas(BaseFold,RelPath,DirectFold,
                  canvas_width=500, canvas_height=500):
        '''调整Canvas Size，并放缩图片，让某一边达到canvas长度'''
        for filename in RelPath:
            old_image_path = BaseFold+filename
            new_image_path = DirectFold+filename
            im = Image.open(old_image_path)
            old_width, old_height = im.size

            # Rescale Pic
            hsize = old_height/old_width
            if canvas_width/old_width > canvas_height/old_height:
                basewidth = canvas_width
            else:
                basewidth = canvas_height/hsize
            im = im.resize((int(basewidth),int(basewidth*hsize)),Image.ANTIALIAS)
            old_width, old_height = im.size

            # Center the image
            x1 = int(math.floor((canvas_width - old_width) / 2))
            y1 = int(math.floor((canvas_height - old_height) / 2))

            mode = im.mode
            if len(mode) == 1:  # L, 1
                new_background = (255)
            if len(mode) == 3:  # RGB
                new_background = (255, 255, 255)
            if len(mode) == 4:  # RGBA, CMYK
                new_background = (255, 255, 255, 255)
            
            newImage = Image.new(mode, (canvas_width, canvas_height), new_background)
            newImage.paste(im, (x1, y1, x1 + old_width, y1 + old_height))
            newImage.save(new_image_path)
            im.close()
            newImage.close()

    def Produce_Model(BaseFold,RelPath,canvas_width=500, canvas_height=500,TagPath = '',Tag=1):
        '''以TagPath为模板生成类似的图片'''
        if TagPath:
            img_ = Image.open(TagPath)
            x0,x1,y0,y1 = Pic.Get_Boundary(img_)
            canvas_width, canvas_height = img_.size
            mode = img_.mode
            if len(mode) == 1:  # L, 1
                new_background = (255)
            if len(mode) == 3:  # RGB
                new_background = (255, 255, 255)
            if len(mode) == 4:  # RGBA, CMYK
                new_background = (255, 255, 255, 255)
            img_.close()

            for filename in RelPath:
                old_image_path = BaseFold+filename
                new_image_path = BaseFold+'New-'+filename
                im = Image.open(old_image_path)
                old_width, old_height = im.size

                # Rescale Pic
                tag_width, tag_height = y1-y0, x1-x0
                hsize = old_height/old_width
                ## 按照width来放缩
                if Tag:
                    basewidth = tag_width*Tag
                else:
                    basewidth = tag_height/hsize
                im = im.resize((int(basewidth),int(basewidth*hsize)),Image.ANTIALIAS)
                old_width, old_height = im.size

                ## 新建同样大小的画布，粘贴上scale过的图片
                newImage = Image.new(mode, (canvas_width, canvas_height), new_background)
                x0 = int((canvas_height-old_height)/2)
                y0 = int((canvas_width-old_width)/2)
                newImage.paste(im, (y0, x0, y0+old_width, x0+old_height))
                newImage.save(new_image_path)
                im.close()
                newImage.close()

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