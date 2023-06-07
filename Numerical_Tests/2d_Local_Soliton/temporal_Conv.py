import sys
import os
os.chdir('/Users/hjs/Documents/2022/Collocation_NLS/Draw_Pic')
sys.path.append('/Users/hjs/Documents/2022/Collocation_NLS/Draw_Pic')

import matplotlib.pyplot as plt
import Data
from Package_MyCode import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle()

# Read Convergece Errors and Plot
Res = Data.LocalSoliton_Temporal
Legendlabel = Res['Params']['Legendlabel']


_order_key,_add_order = Data.GetOrder(Res)
Zoom_xy_scale = [
    {'x': [1/2,2], 'y': [1/3,3]},
    {'x': [1/3,3], 'y': [1/4,4]},
    {'x': [1/4,4], 'y': [1/3,3]}
]
Val_List = [
    [[4e-2,6e-2,1e-1], [1e-3,1e-2]],
    [[5e-2,0.1,5e-1],  [1e-5,1e-4,1e-3]],
    [[1e-2,1e-1,1],    [1e-5,1e-4,1e-3]]
]
XYTick_scale = Data.GenerateXYTick(Val_List)
Zoom=dict(zip(_order_key,Zoom_xy_scale))
XYTick=dict(zip(_order_key,XYTick_scale))

proposed_md = '{}={}'.format(Legendlabel,Res['Params'][Legendlabel])
Legend_scale = [
    [proposed_md,'$\mathcal{O}(\\tau^3)$'],
    [proposed_md,'$\mathcal{O}(\\tau^4)$'],
    [proposed_md,'$\mathcal{O}(\\tau^5)$'],
]
Legend=dict(zip(_order_key,Legend_scale))

ERR_Plot = CanonicalErrPlot(Res, adj_p=[0.7,0.75,0.85], Legend=Legend, Zoom=Zoom)

ERR_Plot.ReadParam()
fig = plt.figure(figsize=(5*ERR_Plot.n_order,5))
ERR_Plot.DrawParam(fig)
# plt.show()

SavePath = '/Users/hjs/Documents/2022/Collocation_NLS/Draw_Pic/Pic'
fig.savefig(os.path.join(SavePath,'{}.eps'.format(Res['Params']['Name'])),dpi=600,format='eps')