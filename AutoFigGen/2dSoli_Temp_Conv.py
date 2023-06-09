import json
import os
import matplotlib.pyplot as plt
from IO_Func import *

MyplotObj =  PlotLineStyleObj()
MyplotObj.SetColorCycle()
MyplotObj.SetLinestyleCycle()
MyplotObj.SetMarkerCycle()
MyplotObj.SetPropCycle()

f_path = '../Numerical_Tests/2d_Soliton/Temp_Conv/Sum_Temp_Conv.json'
# Read Convergece Errors and Plot
with open(f_path, "r") as json_file:
    Res = json.load(json_file)
Res['Params'] = {
        'Legendlabel': 'h',
        'xlabel': '$\\tau$',
        'T': 1,
        'h': '$5\\times2^{-6}$',   # 2*L/256
        'L': 1,
        'N_iter_tol': 1e-8,
        'Name': 'LocalSoliton_Temporal',
        'Param': 'param',
        'Err': 'err',
        'Param_Total': 1,  # final time
        'add_order': 1,
        'var': 'k' # for collocation
    }

Legendlabel = Res['Params']['Legendlabel']

_order_key = ["2","3","4"]
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
XYTick_scale = GenerateXYTick(Val_List)
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
plt.show()